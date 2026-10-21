// Rust implementations of coverage.

#[derive(Clone, Debug, Default)]
pub struct CovRun {
    pub seq: i64,
    pub start: i64,
    pub end: i64,
    pub depth: i64,
}

fn cov_key(minus: bool, seq: i64, pos: i64, end: bool) -> u64 {
    if !(0..1 << 30).contains(&seq) {
        panic!("coverage keys: genome sequence index {} is outside 0 .. 2^30", seq);
    }
    if !(0..1 << 32).contains(&pos) {
        panic!("coverage keys: genome position {} is outside 0 .. 2^32", pos);
    }
    (minus as u64) << 63 | (seq as u64) << 33 | (pos as u64) << 1 | end as u64
}

// Keys for the genome image of transcript bases [a, b).
fn cov_add(m: &BioMapping<BioTranscript, BioGenome>, a: i64, b: i64, out: &mut Vec<u64>) {
    for blk in &m.blocks {
        let (s, e) = (blk.srcStart as i64, (blk.srcStart + blk.len) as i64);
        let (lo, hi) = (a.max(s), b.min(e));
        if lo >= hi {
            continue;
        }
        let d = blk.dstStart as i64;
        let (gs, ge, minus) = match blk.orient {
            BioOrient::Fwd => (d + (lo - s), d + (hi - s), false),
            BioOrient::Rev => (d + (e - hi), d + (e - lo), true),
        };
        out.push(cov_key(minus, m.dst, gs, false));
        out.push(cov_key(minus, m.dst, ge, true));
    }
}

fn cov_keys(refs: &BioRefModels, lens: &[(i64, i64)], xs: &[BioFragHits]) -> Vec<u64> {
    if lens.len() != xs.len() {
        panic!("coverage: {} reads but {} placements", lens.len(), xs.len());
    }
    let mut out = Vec::new();
    for (f, &(l1, l2)) in xs.iter().zip(lens) {
        let mut best: Option<&BioTxHit> = None;
        for h in &f.hits {
            if best.map_or(true, |b| h.weight > b.weight) {
                best = Some(h);
            }
        }
        let Some(h) = best else { continue };
        let Some(m) = refs.maps.get(h.tx as usize) else { continue };
        if m.dst < 0 {
            continue;
        }
        match h.mates {
            0 if h.fragLen > 0 => {
                let (left, right) = if h.fwd { (l1, l2) } else { (l2, l1) };
                let end = h.pos + h.fragLen;
                let rs = (end - right).max(h.pos + left);
                cov_add(m, h.pos, h.pos + left, &mut out);
                cov_add(m, rs, end, &mut out);
            }
            2 => cov_add(m, h.pos, h.pos + l2, &mut out),
            _ => cov_add(m, h.pos, h.pos + l1, &mut out),
        }
    }
    out
}

pub fn cov_keys_pairs(_anchor: i64, refs: &BioRefModels, reads: &Vec<(BioFastqEntry, BioFastqEntry)>, hits: &Vec<BioFragHits>) -> Vec<u64> {
    let lens: Vec<(i64, i64)> = reads.iter().map(|(a, b)| (a.bases.bases.len() as i64, b.bases.bases.len() as i64)).collect();
    cov_keys(refs, &lens, hits)
}

pub fn cov_keys_reads(_anchor: i64, refs: &BioRefModels, reads: &Vec<BioFastqEntry>, hits: &Vec<BioFragHits>) -> Vec<u64> {
    let lens: Vec<(i64, i64)> = reads.iter().map(|a| (a.bases.bases.len() as i64, 0)).collect();
    cov_keys(refs, &lens, hits)
}

// bigWig (UCSC bbi format, version 4, no zoom levels)
// ---------------------------------------------------

const BW_ITEMS_PER_SLOT: usize = 1024;
const BW_BLOCK_SIZE: usize = 256;

fn bw_u16(b: &mut Vec<u8>, x: u16) {
    b.extend_from_slice(&x.to_le_bytes());
}
fn bw_u32(b: &mut Vec<u8>, x: u32) {
    b.extend_from_slice(&x.to_le_bytes());
}
fn bw_u64(b: &mut Vec<u8>, x: u64) {
    b.extend_from_slice(&x.to_le_bytes());
}

// One index entry: (start chrom, start, end chrom, end, offset, size).
type BwSpan = (u32, u32, u32, u32, u64, u64);

// The R-tree over the sections, root first; offsets are absolute.
fn bw_index(items: &[BwSpan], at: u64, end_offset: u64) -> Vec<u8> {
    let mut out = Vec::new();
    bw_u32(&mut out, 0x2468ACE0);
    bw_u32(&mut out, BW_BLOCK_SIZE as u32);
    bw_u64(&mut out, items.len() as u64);
    let (first, last) = (items.first().copied().unwrap_or_default(), items.last().copied().unwrap_or_default());
    bw_u32(&mut out, first.0);
    bw_u32(&mut out, first.1);
    bw_u32(&mut out, last.2);
    bw_u32(&mut out, last.3);
    bw_u64(&mut out, end_offset);
    bw_u32(&mut out, BW_ITEMS_PER_SLOT as u32);
    bw_u32(&mut out, 0);
    // Levels bottom up: each node covers up to BW_BLOCK_SIZE children.
    let mut levels: Vec<Vec<BwSpan>> = vec![items.to_vec()];
    while levels.last().unwrap().len() > BW_BLOCK_SIZE {
        let up: Vec<BwSpan> = levels
            .last()
            .unwrap()
            .chunks(BW_BLOCK_SIZE)
            .map(|c| (c[0].0, c[0].1, c[c.len() - 1].2, c[c.len() - 1].3, 0, 0))
            .collect();
        levels.push(up);
    }
    // Nodes per level, root level first, and their sizes and offsets.
    let node_counts: Vec<usize> = levels.iter().rev().map(|l| (l.len() + BW_BLOCK_SIZE - 1) / BW_BLOCK_SIZE).map(|n| n.max(1)).collect();
    let mut offsets: Vec<Vec<u64>> = Vec::new();
    let mut pos = at + out.len() as u64;
    for (depth, l) in levels.iter().rev().enumerate() {
        let leaf = depth == levels.len() - 1;
        let item = if leaf { 32 } else { 24 };
        let mut offs = Vec::new();
        for c in 0..node_counts[depth] {
            offs.push(pos);
            let n = l.len().saturating_sub(c * BW_BLOCK_SIZE).min(BW_BLOCK_SIZE);
            pos += 4 + (n * item) as u64;
        }
        offsets.push(offs);
    }
    for (depth, l) in levels.iter().rev().enumerate() {
        let leaf = depth == levels.len() - 1;
        let chunks: Vec<&[BwSpan]> = if l.is_empty() { vec![&[][..]] } else { l.chunks(BW_BLOCK_SIZE).collect() };
        for (c, chunk) in chunks.iter().enumerate() {
            out.push(leaf as u8);
            out.push(0);
            bw_u16(&mut out, chunk.len() as u16);
            for (k, s) in chunk.iter().enumerate() {
                bw_u32(&mut out, s.0);
                bw_u32(&mut out, s.1);
                bw_u32(&mut out, s.2);
                bw_u32(&mut out, s.3);
                if leaf {
                    bw_u64(&mut out, s.4);
                    bw_u64(&mut out, s.5);
                } else {
                    bw_u64(&mut out, offsets[depth + 1][c * BW_BLOCK_SIZE + k]);
                }
            }
        }
    }
    out
}

// The chromosome B+ tree: keys sorted by name, ids are indices into `names`.
fn bw_chrom_tree(names: &[String], lens: &[u32], at: u64) -> Vec<u8> {
    let key = names.iter().map(|s| s.len()).max().unwrap_or(1).max(1);
    let mut items: Vec<(Vec<u8>, u32, u32)> = names
        .iter()
        .enumerate()
        .map(|(k, n)| {
            let mut b = n.as_bytes().to_vec();
            b.resize(key, 0);
            (b, k as u32, lens[k])
        })
        .collect();
    items.sort_by(|a, b| a.0.cmp(&b.0));
    if items.windows(2).any(|w| w[0].0 == w[1].0) {
        panic!("bigWig: sequence names repeat");
    }
    let block = items.len().clamp(1, BW_BLOCK_SIZE);
    let mut out = Vec::new();
    bw_u32(&mut out, 0x78CA8C91);
    bw_u32(&mut out, block as u32);
    bw_u32(&mut out, key as u32);
    bw_u32(&mut out, 8);
    bw_u64(&mut out, items.len() as u64);
    bw_u64(&mut out, 0);
    // Levels bottom up: the first key of each node of the level below.
    let mut levels: Vec<Vec<Vec<u8>>> = vec![items.iter().map(|x| x.0.clone()).collect()];
    while levels.last().unwrap().len() > block {
        let up: Vec<Vec<u8>> = levels.last().unwrap().chunks(block).map(|c| c[0].clone()).collect();
        levels.push(up);
    }
    let node_counts: Vec<usize> = levels.iter().rev().map(|l| ((l.len() + block - 1) / block).max(1)).collect();
    let mut offsets: Vec<Vec<u64>> = Vec::new();
    let mut pos = at + out.len() as u64;
    for (depth, l) in levels.iter().rev().enumerate() {
        let leaf = depth == levels.len() - 1;
        let item = key + 8;
        let mut offs = Vec::new();
        for c in 0..node_counts[depth] {
            offs.push(pos);
            let n = l.len().saturating_sub(c * block).min(block);
            pos += 4 + (n * item) as u64;
        }
        offsets.push(offs);
    }
    for (depth, l) in levels.iter().rev().enumerate() {
        let leaf = depth == levels.len() - 1;
        let chunks: Vec<&[Vec<u8>]> = if l.is_empty() { vec![&[][..]] } else { l.chunks(block).collect() };
        for (c, chunk) in chunks.iter().enumerate() {
            out.push(leaf as u8);
            out.push(0);
            bw_u16(&mut out, chunk.len() as u16);
            for (k, kb) in chunk.iter().enumerate() {
                out.extend_from_slice(kb);
                if leaf {
                    let it = &items[c * block + k];
                    bw_u32(&mut out, it.1);
                    bw_u32(&mut out, it.2);
                } else {
                    bw_u64(&mut out, offsets[depth + 1][c * block + k]);
                }
            }
        }
    }
    out
}

// A bigWig written as its runs arrive, in order of sequence index and
// position; only the section index is kept.
struct BwWriter {
    path: String,
    file: std::io::BufWriter<std::fs::File>,
    at: u64,
    chrom_offset: u64,
    data_offset: u64,
    summary_offset: u64,
    lens: Vec<u32>,
    sec: Vec<(u32, u32, u32, f32)>,
    index: Vec<BwSpan>,
    max_raw: usize,
    stats: (u64, f64, f64, f64, f64),
}

impl BwWriter {
    fn put(&mut self, b: &[u8]) {
        use std::io::Write;
        let path = &self.path;
        self.file.write_all(b).unwrap_or_else(|e| panic!("writeCoverage: cannot write {}: {}", path, e));
        self.at += b.len() as u64;
    }

    fn create(path: &str, names: &[String], lens: &[u32]) -> BwWriter {
        let f = std::fs::File::create(path).unwrap_or_else(|e| panic!("writeCoverage: cannot create {}: {}", path, e));
        let mut w = BwWriter {
            path: path.to_string(),
            file: std::io::BufWriter::with_capacity(1 << 20, f),
            at: 0,
            chrom_offset: 0,
            data_offset: 0,
            summary_offset: 0,
            lens: lens.to_vec(),
            sec: Vec::new(),
            index: Vec::new(),
            max_raw: 0,
            stats: (0, f64::INFINITY, f64::NEG_INFINITY, 0.0, 0.0),
        };
        w.put(&[0u8; 64]);
        w.summary_offset = w.at;
        w.put(&[0u8; 40]);
        w.chrom_offset = w.at;
        let tree = bw_chrom_tree(names, lens, w.at);
        w.put(&tree);
        w.data_offset = w.at;
        w.put(&[0u8; 8]);
        w
    }

    fn add(&mut self, seq: u32, s: u32, e: u32, v: f32) {
        if e > self.lens[seq as usize] {
            panic!("writeCoverage: coverage past the end of sequence {} in {}", seq, self.path);
        }
        if !self.sec.is_empty() && (self.sec[0].0 != seq || self.sec.len() == BW_ITEMS_PER_SLOT) {
            self.flush_section();
        }
        let w = (e - s) as f64;
        let x = v as f64;
        let t = &mut self.stats;
        t.0 += (e - s) as u64;
        t.1 = t.1.min(x);
        t.2 = t.2.max(x);
        t.3 += x * w;
        t.4 += x * x * w;
        self.sec.push((seq, s, e, v));
    }

    fn flush_section(&mut self) {
        use flate2::write::ZlibEncoder;
        use std::io::Write;
        let sec = std::mem::take(&mut self.sec);
        if sec.is_empty() {
            return;
        }
        let mut raw = Vec::with_capacity(24 + 12 * sec.len());
        bw_u32(&mut raw, sec[0].0);
        bw_u32(&mut raw, sec[0].1);
        bw_u32(&mut raw, sec[sec.len() - 1].2);
        bw_u32(&mut raw, 0);
        bw_u32(&mut raw, 0);
        raw.push(1);
        raw.push(0);
        bw_u16(&mut raw, sec.len() as u16);
        for r in &sec {
            bw_u32(&mut raw, r.1);
            bw_u32(&mut raw, r.2);
            raw.extend_from_slice(&r.3.to_le_bytes());
        }
        self.max_raw = self.max_raw.max(raw.len());
        let mut z = ZlibEncoder::new(Vec::new(), flate2::Compression::default());
        z.write_all(&raw).unwrap();
        let packed = z.finish().unwrap();
        let off = self.at;
        self.put(&packed);
        self.index.push((sec[0].0, sec[0].1, sec[0].0, sec[sec.len() - 1].2, off, packed.len() as u64));
    }

    fn finish(mut self) {
        use std::io::{Seek, Write};
        self.flush_section();
        let index_offset = self.at;
        let idx = bw_index(&self.index, index_offset, index_offset);
        self.put(&idx);
        let path = self.path.clone();
        let fail = |e: std::io::Error| -> ! { panic!("writeCoverage: cannot write {}: {}", path, e) };
        let mut f = self.file.into_inner().unwrap_or_else(|e| fail(e.into_error()));
        let mut h = Vec::with_capacity(64);
        bw_u32(&mut h, 0x888FFC26);
        bw_u16(&mut h, 4);
        bw_u16(&mut h, 0);
        bw_u64(&mut h, self.chrom_offset);
        bw_u64(&mut h, self.data_offset);
        bw_u64(&mut h, index_offset);
        bw_u16(&mut h, 0);
        bw_u16(&mut h, 0);
        bw_u64(&mut h, 0);
        bw_u64(&mut h, self.summary_offset);
        bw_u32(&mut h, self.max_raw.max(1) as u32);
        bw_u64(&mut h, 0);
        let (n, mut mn, mut mx, sum, sq) = self.stats;
        if n == 0 {
            mn = 0.0;
            mx = 0.0;
        }
        let mut summary = Vec::with_capacity(40);
        bw_u64(&mut summary, n);
        summary.extend_from_slice(&mn.to_le_bytes());
        summary.extend_from_slice(&mx.to_le_bytes());
        summary.extend_from_slice(&sum.to_le_bytes());
        summary.extend_from_slice(&sq.to_le_bytes());
        let mut count = Vec::with_capacity(8);
        bw_u64(&mut count, self.index.len() as u64);
        for (off, b) in [(0u64, &h), (self.summary_offset, &summary), (self.data_offset, &count)] {
            f.seek(std::io::SeekFrom::Start(off)).unwrap_or_else(|e| fail(e));
            f.write_all(b).unwrap_or_else(|e| fail(e));
        }
        f.sync_all().unwrap_or_else(|e| fail(e));
    }
}

// Maximal runs of equal nonzero depth from depth changes ordered by strand,
// sequence and position: changes at one position are summed before the
// depth moves, and a run is emitted when the depth changes.
struct CovRunBuilder {
    cur: Option<(usize, i64)>,
    at: i64,
    net: i64,
    depth: i64,
    run: Option<i64>,
    out: [Vec<CovRun>; 2],
}

impl CovRunBuilder {
    fn new() -> CovRunBuilder {
        CovRunBuilder { cur: None, at: 0, net: 0, depth: 0, run: None, out: [Vec::new(), Vec::new()] }
    }

    fn close(&mut self) {
        let Some((strand, seq)) = self.cur else { return };
        let next = self.depth + self.net;
        if next < 0 {
            panic!("coverage on sequence {} falls below 0", seq);
        }
        if next != self.depth {
            if let Some(a) = self.run.take() {
                self.out[strand].push(CovRun { seq, start: a, end: self.at, depth: self.depth });
            }
            if next > 0 {
                self.run = Some(self.at);
            }
            self.depth = next;
        }
    }

    fn add(&mut self, strand: usize, seq: i64, pos: i64, delta: i64) {
        if self.cur == Some((strand, seq)) && pos == self.at {
            self.net += delta;
            return;
        }
        if let Some(c) = self.cur {
            if (strand, seq, pos) < (c.0, c.1, self.at) {
                panic!("coverage changes are not in order");
            }
        }
        self.close();
        if self.cur != Some((strand, seq)) {
            if self.depth != 0 {
                panic!("coverage on sequence {} does not return to 0", self.cur.unwrap().1);
            }
            self.cur = Some((strand, seq));
        }
        self.at = pos;
        self.net = delta;
    }

    fn finish(&mut self) {
        self.close();
        if self.depth != 0 {
            panic!("coverage does not return to 0");
        }
    }

    fn take(&mut self, strand: usize) -> Vec<CovRun> {
        std::mem::take(&mut self.out[strand])
    }
}

fn cov_drain<T, P: Fn() -> (bool, String, Vec<T>)>(what: &str, pull: &P) -> Option<Vec<T>> {
    let (ok, msg, xs) = pull();
    if !ok {
        panic!("{}: {}", what, msg);
    }
    if xs.is_empty() { None } else { Some(xs) }
}

pub fn cov_runs<A, B, P>(plus: A, minus: B, pull: P)
where
    A: Fn(&Vec<CovRun>),
    B: Fn(&Vec<CovRun>),
    P: Fn() -> (bool, String, Vec<u64>),
{
    let mut b = CovRunBuilder::new();
    let mut prev = 0u64;
    let flush = |b: &mut CovRunBuilder| {
        let p = b.take(0);
        if !p.is_empty() {
            plus(&p);
        }
        let m = b.take(1);
        if !m.is_empty() {
            minus(&m);
        }
    };
    while let Some(xs) = cov_drain("coverageRuns", &pull) {
        for &k in &xs {
            if k < prev {
                panic!("coverageRuns: keys are not ascending");
            }
            prev = k;
            let strand = (k >> 63) as usize;
            let seq = ((k >> 33) & ((1 << 30) - 1)) as i64;
            let pos = ((k >> 1) & 0xFFFF_FFFF) as i64;
            b.add(strand, seq, pos, if k & 1 == 1 { -1 } else { 1 });
        }
        flush(&mut b);
    }
    b.finish();
    flush(&mut b);
}

// The lengths of `names`, looked up in (ref_names, ref_lens).
fn cov_lens(names: &[String], ref_names: &[String], ref_lens: &[i64]) -> Vec<u32> {
    let len_of: std::collections::HashMap<&str, i64> = ref_names.iter().map(|n| n.as_str()).zip(ref_lens.iter().copied()).collect();
    names
        .iter()
        .map(|n| {
            let l = *len_of.get(n.as_str()).unwrap_or_else(|| panic!("coverage tracks: no length for sequence '{}'", n));
            u32::try_from(l).unwrap_or_else(|_| panic!("coverage tracks: sequence '{}' has length {}, outside a bigWig's 0 .. 2^32", n, l))
        })
        .collect()
}

struct CovTrack<'a> {
    names: &'a [String],
    bed: std::io::BufWriter<std::fs::File>,
    bed_path: String,
    bw: BwWriter,
}

impl<'a> CovTrack<'a> {
    fn create(names: &'a [String], lens: &[u32], bed: &str, bw: &str) -> CovTrack<'a> {
        CovTrack {
            names,
            bed: std::io::BufWriter::with_capacity(1 << 20, std::fs::File::create(bed).unwrap_or_else(|e| panic!("coverage tracks: cannot create {}: {}", bed, e))),
            bed_path: bed.to_string(),
            bw: BwWriter::create(bw, names, lens),
        }
    }

    fn add(&mut self, r: &CovRun) {
        use std::io::Write;
        let name = self.names.get(r.seq as usize).unwrap_or_else(|| panic!("coverage tracks: no name for sequence {}", r.seq));
        writeln!(self.bed, "{}\t{}\t{}\t{}", name, r.start, r.end, r.depth).unwrap_or_else(|e| panic!("coverage tracks: cannot write {}: {}", self.bed_path, e));
        self.bw.add(r.seq as u32, r.start as u32, r.end as u32, r.depth as f32);
    }

    fn finish(self) {
        use std::io::Write;
        let mut bed = self.bed;
        bed.flush().unwrap_or_else(|e| panic!("coverage tracks: cannot write {}: {}", self.bed_path, e));
        self.bw.finish();
    }
}

pub fn cov_write_track<P>(names: &Vec<String>, ref_names: &Vec<String>, ref_lens: &Vec<i64>, bed: &String, bw: &String, pull: P)
where
    P: Fn() -> (bool, String, Vec<CovRun>),
{
    let lens = cov_lens(names, ref_names, ref_lens);
    let mut t = CovTrack::create(names, &lens, bed, bw);
    while let Some(xs) = cov_drain("writeCoverageTrack", &pull) {
        for r in &xs {
            t.add(r);
        }
    }
    t.finish();
}

// The runs of one strand as depth changes, a batch pulled at a time.
struct CovChanges<P> {
    pull: P,
    buf: std::collections::VecDeque<(i64, i64, i64)>,
    done: bool,
}

impl<P: Fn() -> (bool, String, Vec<CovRun>)> CovChanges<P> {
    fn peek(&mut self) -> Option<(i64, i64, i64)> {
        while self.buf.is_empty() && !self.done {
            match cov_drain("writeCombinedTrack", &self.pull) {
                None => self.done = true,
                Some(xs) => {
                    for r in xs {
                        self.buf.push_back((r.seq, r.start, r.depth));
                        self.buf.push_back((r.seq, r.end, -r.depth));
                    }
                }
            }
        }
        self.buf.front().copied()
    }
}

pub fn cov_write_combined<P, Q>(names: &Vec<String>, ref_names: &Vec<String>, ref_lens: &Vec<i64>, bed: &String, bw: &String, plus: P, minus: Q)
where
    P: Fn() -> (bool, String, Vec<CovRun>),
    Q: Fn() -> (bool, String, Vec<CovRun>),
{
    let lens = cov_lens(names, ref_names, ref_lens);
    let mut t = CovTrack::create(names, &lens, bed, bw);
    let mut a = CovChanges { pull: plus, buf: Default::default(), done: false };
    let mut c = CovChanges { pull: minus, buf: Default::default(), done: false };
    let mut b = CovRunBuilder::new();
    loop {
        let (x, y) = (a.peek(), c.peek());
        let take_a = match (x, y) {
            (None, None) => break,
            (Some(_), None) => true,
            (None, Some(_)) => false,
            (Some(p), Some(q)) => (p.0, p.1) <= (q.0, q.1),
        };
        let (seq, pos, delta) = if take_a { a.buf.pop_front().unwrap() } else { c.buf.pop_front().unwrap() };
        b.add(0, seq, pos, delta);
        for r in b.take(0) {
            t.add(&r);
        }
    }
    b.finish();
    for r in b.take(0) {
        t.add(&r);
    }
    t.finish();
}
