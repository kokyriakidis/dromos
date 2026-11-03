//! The gbz crate behind a C interface, for dromos: the crate reads the GBZ (its headers and tags, the path metadata,
//! GBWT path walks, node sequences, edges and the bidirectional search) and dromos does everything above it in Mojo
//! (dromos/graph.mojo wraps this). A thin layer: no pgphase or pgcatalog logic lives here.
//!
//! Conventions:
//! - A call that can fail returns 0, or -1 with a message `dgbz_error` reads (per thread, until its next failure).
//! - Strings go in as a pointer and a length, UTF-8, no terminator.
//! - A node visit is a u64 as the GBWT encodes it: `node << 1 | reverse`.
//! - Results of unknown length come back in a `DgbzU64s` or `DgbzBytes` the caller frees with its `_free`; borrowed
//!   pointers (node sequences, names) stay valid as long as the handle they come from.
//! - Handles are read-only once made: any thread may call into one at the same time as others.

use gbz::headers::{GBWTPayload, GBZPayload, Header, SequencesPayload};
use gbz::support::{Orientation, StringArray, Tags};
use gbz::bwt::Record;
use gbz::{BidirectionalState, Metadata, SearchState, ENDMARKER, GBWT, GBZ};
use simple_sds::int_vector::IntVector;
use simple_sds::ops::{BitVec, Select, Vector};
use simple_sds::raw_vector::{AccessRaw, RawVector};
use simple_sds::serialize::Serialize;
use simple_sds::sparse_vector::SparseVector;
use std::cell::RefCell;
use std::fs::File;
use std::io::{BufReader, Seek, SeekFrom};
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::sync::OnceLock;

/// The interface's version: dromos checks it, so a library left from another build is not used by mistake.
pub const ABI: u32 = 2;

#[no_mangle]
pub extern "C" fn dgbz_abi() -> u32 {
    ABI
}

thread_local! {
    static ERROR: RefCell<Vec<u8>> = const { RefCell::new(Vec::new()) };
}

fn fail(message: String) -> i32 {
    ERROR.with(|e| *e.borrow_mut() = message.into_bytes());
    -1
}

/// The last failure's message on this thread (its length into `len`).
#[no_mangle]
pub unsafe extern "C" fn dgbz_error(len: *mut usize) -> *const u8 {
    ERROR.with(|e| {
        let e = e.borrow();
        *len = e.len();
        e.as_ptr()
    })
}

unsafe fn text<'a>(ptr: *const u8, len: usize) -> Result<&'a str, String> {
    if len == 0 {
        return Ok("");
    }
    std::str::from_utf8(std::slice::from_raw_parts(ptr, len)).map_err(|_| String::from("a name that is not UTF-8"))
}

/// Runs `f`, a panic in the crate becoming a failure.
fn guarded(f: impl FnOnce() -> Result<(), String>) -> i32 {
    match catch_unwind(AssertUnwindSafe(f)) {
        Ok(Ok(())) => 0,
        Ok(Err(message)) => fail(message),
        Err(_) => fail(String::from("the gbz crate panicked")),
    }
}

fn orientation(reverse: i32) -> Orientation {
    if reverse != 0 {
        Orientation::Reverse
    } else {
        Orientation::Forward
    }
}

fn encode(node: usize, o: Orientation) -> u64 {
    gbz::support::encode_node(node, o) as u64
}

// Buffers handed to the caller.

#[repr(C)]
pub struct DgbzU64s {
    pub ptr: *mut u64,
    pub len: usize,
}

#[repr(C)]
pub struct DgbzBytes {
    pub ptr: *mut u8,
    pub len: usize,
}

unsafe fn give<T>(v: Vec<T>) -> (*mut T, usize) {
    let b = v.into_boxed_slice();
    let len = b.len();
    (Box::into_raw(b) as *mut T, len)
}

unsafe fn take<T>(ptr: *mut T, len: usize) {
    if !ptr.is_null() {
        drop(Box::from_raw(std::ptr::slice_from_raw_parts_mut(ptr, len)));
    }
}

#[no_mangle]
pub unsafe extern "C" fn dgbz_u64s_free(v: *mut DgbzU64s) {
    take((*v).ptr, (*v).len);
    (*v).ptr = std::ptr::null_mut();
    (*v).len = 0;
}

#[no_mangle]
pub unsafe extern "C" fn dgbz_bytes_free(v: *mut DgbzBytes) {
    take((*v).ptr, (*v).len);
    (*v).ptr = std::ptr::null_mut();
    (*v).len = 0;
}

// Reads from the file's start, without loading the graph (the GBZ's own pieces, read in its order).

fn open(path: &str) -> Result<BufReader<File>, String> {
    let f = File::open(path).map_err(|e| format!("{path}: {e}"))?;
    Ok(BufReader::with_capacity(1 << 20, f))
}

/// Reads the GBZ up to the end of the GBWT's tags; returns them.
fn gbwt_tags(r: &mut BufReader<File>, path: &str) -> Result<Tags, String> {
    let err = |what: &str, e: std::io::Error| format!("{path}: {what}: {e}");
    Header::<GBZPayload>::load(r).map_err(|e| err("not a GBZ file", e))?;
    Tags::load(r).map_err(|e| err("GBZ tags", e))?;
    Header::<GBWTPayload>::load(r).map_err(|e| err("GBWT header", e))?;
    Tags::load(r).map_err(|e| err("GBWT tags", e))
}

/// The GBWT's tag `key` (`found` 0 when it has none), read from the start of the file.
#[no_mangle]
pub unsafe extern "C" fn dgbz_file_gbwt_tag(
    path: *const u8,
    path_len: usize,
    key: *const u8,
    key_len: usize,
    out: *mut DgbzBytes,
    found: *mut i32,
) -> i32 {
    guarded(|| {
        let path = text(path, path_len)?;
        let key = text(key, key_len)?;
        let tags = gbwt_tags(&mut open(path)?, path)?;
        *found = 0;
        if let Some(v) = tags.get(key) {
            let (ptr, len) = give(v.as_bytes().to_vec());
            *out = DgbzBytes { ptr, len };
            *found = 1;
        }
        Ok(())
    })
}

/// The GBWT's path metadata, the GBWT body before it (the BWT index, then its records and the document array
/// samples, two size-prefixed arrays) read past rather than loaded: a null `out` when the GBZ has none. Free it with
/// `dgbz_metadata_free`.
#[no_mangle]
pub unsafe extern "C" fn dgbz_file_metadata(path: *const u8, path_len: usize, out: *mut *mut Metadata) -> i32 {
    guarded(|| {
        let path = text(path, path_len)?;
        let mut r = open(path)?;
        gbwt_tags(&mut r, path)?;
        let err = |what: &str, e: std::io::Error| format!("{path}: {what}: {e}");
        SparseVector::load(&mut r).map_err(|e| err("BWT index", e))?;
        let bytes = usize::load(&mut r).map_err(|e| err("BWT", e))?;
        r.seek(SeekFrom::Current(bytes.div_ceil(8) as i64 * 8)).map_err(|e| err("BWT", e))?;
        let words = usize::load(&mut r).map_err(|e| err("DA samples", e))?;
        r.seek(SeekFrom::Current(words as i64 * 8)).map_err(|e| err("DA samples", e))?;
        let md = Option::<Metadata>::load(&mut r).map_err(|e| err("metadata", e))?;
        *out = md.map_or(std::ptr::null_mut(), |m| Box::into_raw(Box::new(m)));
        Ok(())
    })
}

#[no_mangle]
pub unsafe extern "C" fn dgbz_metadata_free(md: *mut Metadata) {
    if !md.is_null() {
        drop(Box::from_raw(md));
    }
}

// Path metadata, of a graph (`dgbz_metadata`) or read from a file (`dgbz_file_metadata`).

/// A path's name: its sample, contig, phase (haplotype) and fragment (for a reference fragment, its first base).
#[repr(C)]
pub struct DgbzPathName {
    pub sample: u64,
    pub contig: u64,
    pub phase: u64,
    pub fragment: u64,
}

#[no_mangle]
pub unsafe extern "C" fn dgbz_md_paths(md: *const Metadata) -> usize {
    (*md).paths()
}

#[no_mangle]
pub unsafe extern "C" fn dgbz_md_samples(md: *const Metadata) -> usize {
    (*md).samples()
}

#[no_mangle]
pub unsafe extern "C" fn dgbz_md_contigs(md: *const Metadata) -> usize {
    (*md).contigs()
}

/// The names of paths `first .. first + count` into `out`; how many were written (fewer past the last path).
#[no_mangle]
pub unsafe extern "C" fn dgbz_md_path_names(
    md: *const Metadata,
    first: usize,
    count: usize,
    out: *mut DgbzPathName,
) -> usize {
    let md = &*md;
    let mut written = 0;
    for id in first..first.saturating_add(count).min(md.paths()) {
        let Some(p) = md.path(id) else { break };
        *out.add(written) = DgbzPathName {
            sample: p.sample() as u64,
            contig: p.contig() as u64,
            phase: p.phase() as u64,
            fragment: p.fragment as u64,
        };
        written += 1;
    }
    written
}

unsafe fn name(found: Option<&str>, len: *mut usize) -> *const u8 {
    match found {
        Some(s) => {
            *len = s.len();
            s.as_ptr()
        }
        None => {
            *len = 0;
            std::ptr::null()
        }
    }
}

/// Sample `id`'s name, borrowed; null when it has none (the crate's `sample_name` then names it by its number).
#[no_mangle]
pub unsafe extern "C" fn dgbz_md_sample(md: *const Metadata, id: usize, len: *mut usize) -> *const u8 {
    name((*md).sample(id), len)
}

/// Contig `id`'s name, borrowed; null when it has none.
#[no_mangle]
pub unsafe extern "C" fn dgbz_md_contig(md: *const Metadata, id: usize, len: *mut usize) -> *const u8 {
    name((*md).contig(id), len)
}

/// The sample named so; -1 for none.
#[no_mangle]
pub unsafe extern "C" fn dgbz_md_sample_id(md: *const Metadata, name: *const u8, len: usize) -> i64 {
    match text(name, len) {
        Ok(s) => (*md).sample_id(s).map_or(-1, |id| id as i64),
        Err(_) => -1,
    }
}

/// The contig named so; -1 for none.
#[no_mangle]
pub unsafe extern "C" fn dgbz_md_contig_id(md: *const Metadata, name: *const u8, len: usize) -> i64 {
    match text(name, len) {
        Ok(s) => (*md).contig_id(s).map_or(-1, |id| id as i64),
        Err(_) => -1,
    }
}

// A graph, loaded whole by the crate.

/// Loads the GBZ at `path` into `out` (free it with `dgbz_free`). The message of a failure is the crate's.
#[no_mangle]
pub unsafe extern "C" fn dgbz_load(path: *const u8, path_len: usize, out: *mut *mut GBZ) -> i32 {
    guarded(|| {
        let path = text(path, path_len)?;
        let g: GBZ = simple_sds::serialize::load_from(path).map_err(|e| e.to_string())?;
        *out = Box::into_raw(Box::new(g));
        Ok(())
    })
}

#[no_mangle]
pub unsafe extern "C" fn dgbz_free(g: *mut GBZ) {
    if !g.is_null() {
        drop(Box::from_raw(g));
    }
}

/// The graph's path metadata, borrowed; null when it has none.
#[no_mangle]
pub unsafe extern "C" fn dgbz_metadata(g: *const GBZ) -> *const Metadata {
    (*g).metadata().map_or(std::ptr::null(), |m| m as *const Metadata)
}

#[no_mangle]
pub unsafe extern "C" fn dgbz_nodes(g: *const GBZ) -> usize {
    (*g).nodes()
}

#[no_mangle]
pub unsafe extern "C" fn dgbz_min_node(g: *const GBZ) -> usize {
    (*g).min_node()
}

#[no_mangle]
pub unsafe extern "C" fn dgbz_max_node(g: *const GBZ) -> usize {
    (*g).max_node()
}

#[no_mangle]
pub unsafe extern "C" fn dgbz_has_node(g: *const GBZ, node: usize) -> i32 {
    (*g).has_node(node) as i32
}

/// Node `node`'s sequence on its forward strand, borrowed (its length into `len`); null for a node the graph lacks.
#[no_mangle]
pub unsafe extern "C" fn dgbz_sequence(g: *const GBZ, node: usize, len: *mut usize) -> *const u8 {
    match (*g).sequence(node) {
        Some(s) => {
            *len = s.len();
            s.as_ptr()
        }
        None => {
            *len = 0;
            std::ptr::null()
        }
    }
}

/// The sequence lengths of nodes `first .. first + count` into `out`, u64::MAX for an id the graph lacks.
#[no_mangle]
pub unsafe extern "C" fn dgbz_sequence_lengths(g: *const GBZ, first: usize, count: usize, out: *mut u64) {
    let g = &*g;
    let out = std::slice::from_raw_parts_mut(out, count);
    for (i, o) in out.iter_mut().enumerate() {
        *o = g.sequence_len(first + i).map_or(u64::MAX, |n| n as u64);
    }
}

/// The sequences of nodes `first .. first + count`, borrowed: each one's address into `address` (0 for an id the
/// graph lacks) and its length into `len`.
#[no_mangle]
pub unsafe extern "C" fn dgbz_sequences(g: *const GBZ, first: usize, count: usize, address: *mut u64, len: *mut u64) {
    let g = &*g;
    for i in 0..count {
        let (a, n) = g.sequence(first + i).map_or((0, 0), |s| (s.as_ptr() as u64, s.len() as u64));
        *address.add(i) = a;
        *len.add(i) = n;
    }
}

#[no_mangle]
pub unsafe extern "C" fn dgbz_paths(g: *const GBZ) -> usize {
    (*g).paths()
}

/// Path `id`'s visits, walked in the GBWT (`reverse`: its reverse), into `out`.
#[no_mangle]
pub unsafe extern "C" fn dgbz_path(g: *const GBZ, id: usize, reverse: i32, out: *mut DgbzU64s) -> i32 {
    guarded(|| {
        let it = (*g).path(id, orientation(reverse)).ok_or_else(|| format!("no path {id}"))?;
        let visits: Vec<u64> = it.map(|(node, o)| encode(node, o)).collect();
        let (ptr, len) = give(visits);
        *out = DgbzU64s { ptr, len };
        Ok(())
    })
}

unsafe fn edges(it: Option<gbz::gbz::EdgeIter<'_>>, out: *mut u64, cap: usize) -> i64 {
    let Some(it) = it else { return -1 };
    let mut n = 0;
    for (node, o) in it {
        if n < cap {
            *out.add(n) = encode(node, o);
        }
        n += 1;
    }
    n as i64
}

/// The successors of (node, orientation): their count, the first `cap` of them into `out`; -1 for a node the graph
/// lacks.
#[no_mangle]
pub unsafe extern "C" fn dgbz_successors(g: *const GBZ, node: usize, reverse: i32, out: *mut u64, cap: usize) -> i64 {
    edges((*g).successors(node, orientation(reverse)), out, cap)
}

/// The predecessors of (node, orientation), as `dgbz_successors`.
#[no_mangle]
pub unsafe extern "C" fn dgbz_predecessors(
    g: *const GBZ,
    node: usize,
    reverse: i32,
    out: *mut u64,
    cap: usize,
) -> i64 {
    edges((*g).predecessors(node, orientation(reverse)), out, cap)
}

/// A state of bidirectional GBWT search: the last node of each direction (GBWT-encoded) and its range.
#[repr(C)]
#[derive(Clone, Copy)]
pub struct DgbzState {
    pub forward_node: u64,
    pub forward_start: u64,
    pub forward_end: u64,
    pub reverse_node: u64,
    pub reverse_start: u64,
    pub reverse_end: u64,
}

impl DgbzState {
    fn of(s: &BidirectionalState) -> Self {
        DgbzState {
            forward_node: s.forward.node as u64,
            forward_start: s.forward.range.start as u64,
            forward_end: s.forward.range.end as u64,
            reverse_node: s.reverse.node as u64,
            reverse_start: s.reverse.range.start as u64,
            reverse_end: s.reverse.range.end as u64,
        }
    }

    fn state(&self) -> BidirectionalState {
        BidirectionalState {
            forward: SearchState {
                node: self.forward_node as usize,
                range: self.forward_start as usize..self.forward_end as usize,
            },
            reverse: SearchState {
                node: self.reverse_node as usize,
                range: self.reverse_start as usize..self.reverse_end as usize,
            },
        }
    }
}

/// The search state of the haplotypes visiting (node, orientation) into `out`: 1, or 0 when none do.
#[no_mangle]
pub unsafe extern "C" fn dgbz_search_state(g: *const GBZ, node: usize, reverse: i32, out: *mut DgbzState) -> i32 {
    match (*g).search_state(node, orientation(reverse)) {
        Some(s) => {
            *out = DgbzState::of(&s);
            1
        }
        None => 0,
    }
}

unsafe fn states(it: Option<gbz::gbz::StateIter<'_>>, out: *mut DgbzState, cap: usize) -> i64 {
    let Some(it) = it else { return -1 };
    let mut n = 0;
    for s in it {
        if n < cap {
            *out.add(n) = DgbzState::of(&s);
        }
        n += 1;
    }
    n as i64
}

/// The non-empty states `state` extends to forward, one per successor some of its haplotypes take: their count, the
/// first `cap` into `out`; -1 when the state's node is not in the graph.
#[no_mangle]
pub unsafe extern "C" fn dgbz_follow_forward(
    g: *const GBZ,
    state: *const DgbzState,
    out: *mut DgbzState,
    cap: usize,
) -> i64 {
    states((*g).follow_forward(&(*state).state()), out, cap)
}

/// As `dgbz_follow_forward`, extending the state backward.
#[no_mangle]
pub unsafe extern "C" fn dgbz_follow_backward(
    g: *const GBZ,
    state: *const DgbzState,
    out: *mut DgbzState,
    cap: usize,
) -> i64 {
    states((*g).follow_backward(&(*state).state()), out, cap)
}

// A graph's parts, for reading every node and the reference paths fast: the GBWT, loaded by the crate, and the node
// sequences' own pieces, the string array's index, alphabet and packed letters each loaded by the crate (or, in a
// zstd-compressed GBZ, the array itself), read at the same time as the GBWT. The crate's GBZ load unpacks every
// letter on one thread first, most of its time; here letters are unpacked as they are asked for, by node range, on
// the caller's threads. Node ids, real nodes, lengths and letters mean what they mean in the crate's GBZ.

enum Letters {
    /// The string array as stored: string starts, alphabet, letters packed (StringArray::load's pieces), and the
    /// alphabet as a table of 256 for unpacking.
    Packed { starts: SparseVector, alphabet: Vec<u8>, packed: IntVector, table: [u8; 256] },
    /// The string array itself (a zstd-compressed one, StringArray::decompress).
    Strings(StringArray),
}

impl Letters {
    fn count(&self) -> usize {
        match self {
            Letters::Packed { starts, .. } => starts.count_ones(),
            Letters::Strings(s) => s.len(),
        }
    }

    /// The starts of strings `first ..= first + n` (the last may be the end of the letters).
    fn starts(&self, first: usize, n: usize, out: &mut Vec<usize>) {
        out.clear();
        match self {
            Letters::Packed { starts, packed, .. } => {
                out.extend(starts.select_iter(first).take(n + 1).map(|(_, x)| x));
                while out.len() < n + 1 {
                    out.push(packed.len());
                }
            }
            Letters::Strings(s) => {
                let mut at = s.range_len(0..first);
                out.push(at);
                for i in first..first + n {
                    at += s.str_len(i);
                    out.push(at);
                }
            }
        }
    }

    /// The letters of strings `strings` (from letter `start` to `end`) into `out`; false if one is outside the
    /// alphabet.
    fn put(&self, strings: std::ops::Range<usize>, start: usize, end: usize, out: &mut [u8]) -> bool {
        match self {
            Letters::Packed { alphabet, packed, table, .. } => {
                let raw: &RawVector = packed.as_ref();
                let sigma = alphabet.len() as u64;
                let out = &mut out[..end - start];
                match packed.width() {
                    1 => unpack::<1>(raw, start, sigma, table, out),
                    2 => unpack::<2>(raw, start, sigma, table, out),
                    3 => unpack::<3>(raw, start, sigma, table, out),
                    4 => unpack::<4>(raw, start, sigma, table, out),
                    5 => unpack::<5>(raw, start, sigma, table, out),
                    6 => unpack::<6>(raw, start, sigma, table, out),
                    7 => unpack::<7>(raw, start, sigma, table, out),
                    8 => unpack::<8>(raw, start, sigma, table, out),
                    _ => false,
                }
            }
            Letters::Strings(s) => {
                out[..end - start].copy_from_slice(s.range(strings));
                true
            }
        }
    }
}

/// Letters `start ..` of a packed array of `W`-bit letters into `out`, through the alphabet `table`: as
/// IntVector::get reads each (`int(j * W, W)`), a 64-bit window of them per read. False if one is outside the
/// alphabet (`sigma` letters).
fn unpack<const W: usize>(raw: &RawVector, start: usize, sigma: u64, table: &[u8; 256], out: &mut [u8]) -> bool {
    let per = 64 / W;
    let mask = (1u64 << W) - 1;
    let mut over = 0u64;
    let mut j = start;
    let mut chunks = out.chunks_exact_mut(per);
    for chunk in &mut chunks {
        let bits = unsafe { raw.int(j * W, per * W) };
        for (t, o) in chunk.iter_mut().enumerate() {
            let x = (bits >> (t * W)) & mask;
            over |= (x >= sigma) as u64;
            *o = table[x as usize];
        }
        j += per;
    }
    let rest = chunks.into_remainder();
    if !rest.is_empty() {
        let bits = unsafe { raw.int(j * W, rest.len() * W) };
        for (t, o) in rest.iter_mut().enumerate() {
            let x = (bits >> (t * W)) & mask;
            over |= (x >= sigma) as u64;
            *o = table[x as usize];
        }
    }
    over == 0
}

pub struct Parts {
    gbwt: GBWT,
    letters: Letters,
    /// A bit per node sequence: the node is in the graph (the crate's GBZ `real_nodes`), found on first use.
    real: OnceLock<Vec<u64>>,
}

impl Parts {
    /// The sequence of node `node`, if the GBWT's alphabet has it (GBZ::graph_node_to_sequence).
    fn sequence_of(&self, node: usize) -> Option<usize> {
        let g = gbz::support::encode_node(node, Orientation::Forward);
        self.gbwt.has_node(g).then(|| (g - self.gbwt.first_node()) / 2)
    }

    fn min_node(&self) -> usize {
        gbz::support::node_id(self.gbwt.first_node())
    }

    fn max_node(&self) -> usize {
        gbz::support::node_id(self.gbwt.alphabet_size()) - 1
    }

    /// The real nodes, as the crate's GBZ load finds them (the forward nodes whose record is not empty), found on
    /// `threads` threads, by ranges of nodes: a record is empty when `compressed_record` has no edges for it.
    fn real_on(&self, threads: usize) -> &Vec<u64> {
        self.real.get_or_init(|| {
            let n = self.letters.count();
            let mut bits = vec![0u64; n.div_ceil(64)];
            let bwt = self.gbwt.as_ref();
            let first = self.gbwt.first_node();
            // Ranges of 1024 words (65536 nodes), taken in turn by `threads` workers.
            let work = std::sync::Mutex::new(bits.chunks_mut(1024).enumerate().collect::<Vec<_>>());
            std::thread::scope(|s| {
                for _ in 0..threads.max(1) {
                    s.spawn(|| loop {
                        let Some((c, chunk)) = work.lock().unwrap().pop() else { break };
                        for (w, word) in chunk.iter_mut().enumerate() {
                            let s0 = (c * 1024 + w) * 64;
                            let mut x = 0u64;
                            for b in 0..64.min(n.saturating_sub(s0)) {
                                let record = self.gbwt.node_to_record(first + 2 * (s0 + b));
                                if record != ENDMARKER && bwt.compressed_record(record).is_some() {
                                    x |= 1 << b;
                                }
                            }
                            *word = x;
                        }
                    });
                }
            });
            bits
        })
    }

    fn real(&self) -> &Vec<u64> {
        // Once found, without asking how many cores there are (which reads files on Linux).
        match self.real.get() {
            Some(r) => r,
            None => self.real_on(std::thread::available_parallelism().map_or(1, |n| n.get())),
        }
    }

    fn is_real(real: &[u64], s: usize) -> bool {
        real[s / 64] >> (s % 64) & 1 == 1
    }

    /// Whether the graph has node `node` (GBZ::has_node).
    fn has_node(&self, node: usize) -> bool {
        self.sequence_of(node).is_some_and(|s| Parts::is_real(self.real(), s))
    }

    /// The record GBZ::successors reads for (node, orientation); None for a node the graph lacks.
    fn record_of(&self, node: usize, o: Orientation) -> Option<Record<'_>> {
        if !self.has_node(node) {
            return None;
        }
        self.gbwt.as_ref().record(self.gbwt.node_to_record(gbz::support::encode_node(node, o)))
    }

    /// The node range's sequences (first, count) within the GBWT's alphabet.
    fn sequences(&self, first: usize, count: usize) -> (usize, usize, usize) {
        let lo = first.max(self.min_node());
        let hi = (first + count).min(self.max_node() + 1);
        if lo >= hi {
            return (0, 0, 0);
        }
        (lo - first, self.sequence_of(lo).unwrap(), hi - lo)
    }
}

fn load_gbwt(path: &str) -> Result<GBWT, String> {
    let mut r = open(path)?;
    let err = |what: &str, e: std::io::Error| format!("{path}: {what}: {e}");
    let header = Header::<GBZPayload>::load(&mut r).map_err(|e| err("not a GBZ file", e))?;
    header.validate().map_err(|e| format!("{path}: {e}"))?;
    Tags::load(&mut r).map_err(|e| err("GBZ tags", e))?;
    let gbwt = GBWT::load(&mut r).map_err(|e| err("GBWT", e))?;
    if !gbwt.is_bidirectional() {
        return Err(format!("{path}: GBZ: The GBWT index is not bidirectional"));
    }
    Ok(gbwt)
}

fn load_letters(path: &str) -> Result<Letters, String> {
    let mut r = open(path)?;
    gbwt_tags(&mut r, path)?;
    let err = |what: &str, e: std::io::Error| format!("{path}: {what}: {e}");
    SparseVector::load(&mut r).map_err(|e| err("BWT index", e))?;
    let bytes = usize::load(&mut r).map_err(|e| err("BWT", e))?;
    r.seek(SeekFrom::Current(bytes.div_ceil(8) as i64 * 8)).map_err(|e| err("BWT", e))?;
    let words = usize::load(&mut r).map_err(|e| err("DA samples", e))?;
    r.seek(SeekFrom::Current(words as i64 * 8)).map_err(|e| err("DA samples", e))?;
    Option::<Metadata>::load(&mut r).map_err(|e| err("metadata", e))?;
    let header = Header::<SequencesPayload>::load(&mut r).map_err(|e| err("sequences", e))?;
    header.validate().map_err(|e| format!("{path}: {e}"))?;
    if header.version() >= SequencesPayload::ZSTD_VERSION {
        return Ok(Letters::Strings(StringArray::decompress(&mut r).map_err(|e| err("sequences", e))?));
    }
    let starts = SparseVector::load(&mut r).map_err(|e| err("sequences", e))?;
    let alphabet = Vec::<u8>::load(&mut r).map_err(|e| err("sequences", e))?;
    let packed = IntVector::load(&mut r).map_err(|e| err("sequences", e))?;
    if starts.count_ones() > 0 && starts.select(0) != Some(0) {
        return Err(format!("{path}: StringArray: First string does not start at offset 0"));
    }
    if alphabet.len() > 256 {
        return Err(format!("{path}: sequences: an alphabet of more than 256 letters"));
    }
    let mut table = [0u8; 256];
    table[..alphabet.len()].copy_from_slice(&alphabet);
    Ok(Letters::Packed { starts, alphabet, packed, table })
}

/// Loads the graph's parts at `path` into `out` (free it with `dgbz_parts_free`): the GBWT and the node sequences'
/// pieces, read at the same time.
#[no_mangle]
pub unsafe extern "C" fn dgbz_parts_load(path: *const u8, path_len: usize, out: *mut *mut Parts) -> i32 {
    guarded(|| {
        let path = text(path, path_len)?;
        let (gbwt, letters) = std::thread::scope(|s| {
            let letters = s.spawn(|| load_letters(path));
            let gbwt = load_gbwt(path);
            (gbwt, letters.join().unwrap_or_else(|_| Err(String::from("the gbz crate panicked"))))
        });
        let (gbwt, letters) = (gbwt?, letters?);
        if letters.count() != (gbwt.alphabet_size() - gbwt.first_node()) / 2 {
            return Err(format!("{path}: GBZ: Mismatch between GBWT alphabet size and Sequences sequence count"));
        }
        *out = Box::into_raw(Box::new(Parts { gbwt, letters, real: OnceLock::new() }));
        Ok(())
    })
}

#[no_mangle]
pub unsafe extern "C" fn dgbz_parts_free(p: *mut Parts) {
    if !p.is_null() {
        drop(Box::from_raw(p));
    }
}

#[no_mangle]
pub unsafe extern "C" fn dgbz_parts_metadata(p: *const Parts) -> *const Metadata {
    (*p).gbwt.metadata().map_or(std::ptr::null(), |m| m as *const Metadata)
}

#[no_mangle]
pub unsafe extern "C" fn dgbz_parts_max_node(p: *const Parts) -> usize {
    (*p).max_node()
}

/// Finds the real nodes (a scan of the GBWT's records) on `threads` threads, if not found yet: a call that needs
/// them otherwise does it first, on every core.
#[no_mangle]
pub unsafe extern "C" fn dgbz_parts_find_nodes(p: *const Parts, threads: usize) {
    (*p).real_on(threads);
}

#[no_mangle]
pub unsafe extern "C" fn dgbz_parts_paths(p: *const Parts) -> usize {
    (*p).gbwt.sequences() / 2
}

/// Path `id`'s visits, walked in the GBWT (`reverse`: its reverse), into `out`.
#[no_mangle]
pub unsafe extern "C" fn dgbz_parts_path(p: *const Parts, id: usize, reverse: i32, out: *mut DgbzU64s) -> i32 {
    guarded(|| {
        let g = &(*p).gbwt;
        let it = g.sequence(gbz::support::encode_path(id, orientation(reverse)));
        let it = it.ok_or_else(|| format!("no path {id}"))?;
        let (ptr, len) = give(it.map(|x| x as u64).collect::<Vec<u64>>());
        *out = DgbzU64s { ptr, len };
        Ok(())
    })
}

/// As `dgbz_sequence_lengths`.
#[no_mangle]
pub unsafe extern "C" fn dgbz_parts_lengths(p: *const Parts, first: usize, count: usize, out: *mut u64) {
    let p = &*p;
    let real = p.real();
    let out = std::slice::from_raw_parts_mut(out, count);
    out.fill(u64::MAX);
    let (at, s0, n) = p.sequences(first, count);
    let mut starts = Vec::with_capacity(n + 1);
    p.letters.starts(s0, n, &mut starts);
    for i in 0..n {
        if Parts::is_real(real, s0 + i) {
            out[at + i] = (starts[i + 1] - starts[i]) as u64;
        }
    }
}

/// As `dgbz_parts_lengths`, as u32 into `out`: `absent` for an id the graph lacks, a length cut to 32 bits (as
/// `as u32` cuts it). The lengths' total into `sum`; the first id whose length a u32 cannot hold, or -1.
#[no_mangle]
pub unsafe extern "C" fn dgbz_parts_lengths_u32(
    p: *const Parts,
    first: usize,
    count: usize,
    absent: u32,
    out: *mut u32,
    sum: *mut u64,
) -> i64 {
    let p = &*p;
    let real = p.real();
    let out = std::slice::from_raw_parts_mut(out, count);
    out.fill(absent);
    let (at, s0, n) = p.sequences(first, count);
    let mut starts = Vec::with_capacity(n + 1);
    p.letters.starts(s0, n, &mut starts);
    let (mut total, mut long) = (0u64, -1i64);
    for i in 0..n {
        if Parts::is_real(real, s0 + i) {
            let m = (starts[i + 1] - starts[i]) as u64;
            if m > u32::MAX as u64 && long < 0 {
                long = (first + at + i) as i64;
            }
            out[at + i] = m as u32;
            total += m;
        }
    }
    *sum = total;
    long
}

/// The letters of nodes `first .. first + count` that are in the graph, one after another, into `out` (room for
/// `cap`): how many; -1 when they do not fit, -2 for a letter outside the graph's alphabet.
#[no_mangle]
pub unsafe extern "C" fn dgbz_parts_letters(
    p: *const Parts,
    first: usize,
    count: usize,
    out: *mut u8,
    cap: usize,
) -> i64 {
    let p = &*p;
    let real = p.real();
    let out = std::slice::from_raw_parts_mut(out, cap);
    let (_, s0, n) = p.sequences(first, count);
    let mut starts = Vec::with_capacity(n + 1);
    p.letters.starts(s0, n, &mut starts);
    // Runs of nodes in the graph: their letters lie together.
    let mut at = 0;
    let mut i = 0;
    while i < n {
        if !Parts::is_real(real, s0 + i) {
            i += 1;
            continue;
        }
        let mut j = i + 1;
        while j < n && Parts::is_real(real, s0 + j) {
            j += 1;
        }
        let m = starts[j] - starts[i];
        if at + m > cap {
            return -1;
        }
        if !p.letters.put(s0 + i..s0 + j, starts[i], starts[j], &mut out[at..at + m]) {
            return -2;
        }
        at += m;
        i = j;
    }
    at as i64
}

// Edges and bidirectional search on the parts' GBWT, as the crate's GBZ gives them (successors, predecessors,
// search_state, follow_forward), from the same GBWT functions: the same edges and states. Many at once, so a search
// pays one call per batch of states, not one per state.

/// The edge ranks EdgeIter visits in a record: all but an endmarker at rank 0.
fn edge_ranks(record: &Record<'_>) -> std::ops::Range<usize> {
    let first = (record.outdegree() > 0 && record.successor(0) == ENDMARKER) as usize;
    first..record.outdegree()
}

/// For each of the `n` visits at `visits`, its successors (GBZ::successors) or, with `predecessors`, its
/// predecessors (GBZ::predecessors), as visits, one visit's after another into `out` (room for `cap`), and how many
/// into `counts` (-1 for a node the graph lacks). Returns how many visits were done: fewer than `n` when the next
/// one's would not fit.
#[no_mangle]
pub unsafe extern "C" fn dgbz_parts_edges(
    p: *const Parts,
    visits: *const u64,
    n: usize,
    predecessors: i32,
    out: *mut u64,
    cap: usize,
    counts: *mut i64,
) -> usize {
    let p = &*p;
    let mut written = 0;
    for i in 0..n {
        let (node, o) = gbz::support::decode_node(*visits.add(i) as usize);
        // A node's predecessors are the flipped successors of its other orientation.
        let (o, flip) = if predecessors != 0 { (o.flip(), 1) } else { (o, 0) };
        let Some(record) = p.record_of(node, o) else {
            *counts.add(i) = -1;
            continue;
        };
        let ranks = edge_ranks(&record);
        if written + ranks.len() > cap {
            return i;
        }
        for r in ranks.clone() {
            *out.add(written) = record.successor(r) as u64 ^ flip;
            written += 1;
        }
        *counts.add(i) = ranks.len() as i64;
    }
    n
}

/// The search state of the haplotypes visiting (node, orientation) into `out` (GBZ::search_state): 1, or 0 when none
/// do.
#[no_mangle]
pub unsafe extern "C" fn dgbz_parts_search_state(p: *const Parts, node: usize, reverse: i32, out: *mut DgbzState) -> i32 {
    match (*p).gbwt.bd_find(gbz::support::encode_node(node, orientation(reverse))) {
        Some(s) => {
            *out = DgbzState::of(&s);
            1
        }
        None => 0,
    }
}

/// For each of the `n` states at `states`, the non-empty states it extends to forward, one per successor some of
/// its haplotypes take (GBZ::follow_forward), one state's after another into `out` (room for `cap`), and how many
/// into `counts` (-1 when the state's node is not in the graph). Returns how many states were done: fewer than `n`
/// when the next one's might not fit.
#[no_mangle]
pub unsafe extern "C" fn dgbz_parts_follow_forward(
    p: *const Parts,
    states: *const DgbzState,
    n: usize,
    out: *mut DgbzState,
    cap: usize,
    counts: *mut i64,
) -> usize {
    let p = &*p;
    let mut written = 0;
    for i in 0..n {
        let state = (*states.add(i)).state();
        let (node, o) = gbz::support::decode_node(state.forward.node);
        let Some(record) = p.record_of(node, o) else {
            *counts.add(i) = -1;
            continue;
        };
        let ranks = edge_ranks(&record);
        if written + ranks.len() > cap {
            return i;
        }
        let mut k = 0;
        for r in ranks {
            if let Some(s) = GBWT::bd_internal(&record, &state, record.successor(r)) {
                *out.add(written + k) = DgbzState::of(&s);
                k += 1;
            }
        }
        *counts.add(i) = k as i64;
        written += k;
    }
    n
}
