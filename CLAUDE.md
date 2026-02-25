# CLAUDE.md

## Project

Comparing worst-case optimal join (WCOJ / trie-join) with sparse BLAS for
triangle counting. The goal is a pedagogical document showing that the
trie-join framework, when backed by CSR storage, generates code identical
to hand-written CSR triangle counting — and benchmarks confirming this.

## Repository Structure

- `docs/triangle-triejoin-csr.md` — main document (7 sections, iteratively refined)
- `julia/` — Julia benchmark code (standalone, outside RAICode)
- `cpp/` — C++ reference implementation

## Key Context

### The document builds up progressively:
1. Relations as tries (interface, not storage)
2. Two-relation join = SpMM (remarkable that WCOJ subsumes SpMM)
3. Three-relation join: pairwise O(M^2) vs WCOJ O(M^{3/2})
4. Self-join R=S=T with x<y<z
5. Trie-join <-> CSR side by side, with `intersect_neighbors` abstraction
6. In RAICode: TrieState interface -> spans -> CSR (shows generated code = hand-written)
7. Performance results (placeholder tables to be filled)

### Important terminology:
- **CSI** = Contiguous Span Iterator (NOT Column Store Iterator)
- **TSM** = TrieStateMaterialized (B+ tree backed)
- **TSC** = TrieStateContiguous (implicit 1:N keys)
- Trie is an *interface* over data, not a storage format
- Relation = edge list = COO sparse matrix
- Rows are contiguous *within* spans (NOT one span per row)

### Performance results tables:
- **Public table:** RAICode trie, RAICode CSR, julia_csr, cpp_csr
- **Internal table:** incremental improvements: no spans -> +TSM spans -> +TSC spans -> +CSR

### Related RAICode code (in the raicode repo):
- `packages/BerkeleyOperators/src/TrieInterface/trie-state-conjunction.jl`
- `packages/BerkeleyOperators/src/TrieInterface/trie-state-contiguous.jl`
- `bench/Microbench/spans/conjunction-span-bench.jl` — existing triangle benchmark
- `bench/Microbench/spans/spmatvec-bench.jl` — SpMV benchmark (TSC vs TSM)
- `bench/Microbench/spans/SPEC-trie-state-contiguous.md` — TSC specification

## Git

- Remote: `git@github-rai:lums-rai/sparseblas-vs-wcoj-tc.git`
- Uses `github-rai` SSH alias (same as raicode repo)

## Style preferences

- Abstract-first, then concrete: start with trie-join concepts, then show RAICode
- No emojis in documents
- Show side-by-side comparisons (trie-join <-> CSR)
- When showing generated code, annotate each line with the TrieState operation it maps to
