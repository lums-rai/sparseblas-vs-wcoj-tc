# CLAUDE.md

## Project

HPEC 2026 paper: "Generating Fused Graph Kernels from Relational Algebra."
Compares worst-case optimal joins (WCOJ) with sparse linear algebra (GraphBLAS)
for triangle counting. Central thesis: both are declarative frameworks targeting
CSR, but relational algebra generates a *fused* three-way kernel that pairwise
matrix operations cannot express.

## Repository Structure

- `docs/paper/` — HPEC paper (IEEE two-column, Palatino, ~9 pages)
  - `main.tex` — master file, abstract, \input for all sections
  - `refs.bib` — bibliography (~22 entries)
  - `sec-introduction.tex` — Section I: two declarative frameworks framing
  - `sec-background.tex` — Section II: CSR, GraphBLAS (as declarative), AGM bound, variable-at-a-time on CSR
  - `sec-spmm.tex` — Section III: two-relation join = SpMM
  - `sec-triangle.tex` — Section IV: pairwise vs fused, **includes formal WCOJ derivation recipe** (IV-A)
  - `sec-equivalence.tex` — Section V: code equivalence (central result), side-by-side figure, Algorithms 2-3
  - `sec-production.tex` — Section VI: TrieState interface, B+ tree vs CSR backends
  - `sec-lazy.tex` — Section VII: why lazy evaluation can't close the gap (SPORES result)
  - `sec-experiments.tex` — Section VIII: TBD placeholders for benchmarks
  - `sec-related.tex` — Section IX: related work
  - `sec-conclusion.tex` — Section X: conclusion
- `docs/triangle-triejoin-csr.md` — original pedagogical document (reference, not actively edited)
- `julia/` — Julia benchmark code (TBD)
- `cpp/` — C++ reference implementation (TBD)

## Build

```
cd docs/paper && latexmk -pdf main.tex
```

## Key Framing Decisions (established through discussion)

1. **Both frameworks are declarative.** GraphBLAS (sparse LA) and relational
   algebra are both declarative. The paper compares two declarative frameworks,
   not "declarative DB vs imperative HPC."

2. **Loop fusion metaphor.** WCOJ = fused kernel (single pass, simultaneous
   intersection). SpMM = unfused (multiply then filter, two passes). HPC
   audience understands fusion.

3. **Minimal trie formalism.** The trie interface is mentioned in one paragraph
   (Section II-E) as Veldhuizen's formalization, not built up as a standalone
   section. CSR already IS the trie — precomputed index, O(1) row access.

4. **SPORES is more than an aside.** Wang et al. proved RA rewrites are
   complete for LA optimization; LA rewrites are not. This is a dedicated
   section (VII) addressing the "lazy evaluation will fix it" objection.

5. **Formal WCOJ derivation.** Section IV-A gives the three-step recipe for
   constructing a WCOJ evaluation plan from a conjunctive query. Cites
   Veldhuizen (Leapfrog Triejoin), NPRR (Generic Join), FAQ (InsideOut for
   aggregation), and fractional hypertree width.

6. **Variable-at-a-time, not trie ops.** Use N(x), N_R(x) notation instead of
   R.level1, R.open(x). Tables show "Constraints" and "Candidates" columns.

## Important Terminology

- **CSI** = Contiguous Span Iterator (NOT Column Store Iterator)
- **TSM** = TrieStateMaterialized (B+ tree backed)
- **TSC** = TrieStateContiguous (implicit 1:N keys)
- **WCOJ** = Worst-Case Optimal Join
- **FAQ** = Functional Aggregate Queries (Abo Khamis, Ngo, Rudra 2016)
- **AGM bound** = Atserias-Grohe-Marx bound: output ≤ O(M^{3/2}) for triangles
- N(v) = sorted neighbor list of vertex v

## Key References in the Paper

- Veldhuizen 2014 — Leapfrog Triejoin algorithm
- NPRR 2018 — Generic Join, WCOJ theory
- AGM 2013 — output size bounds for conjunctive queries
- Grohe-Marx 2014 — fractional edge covers / hypertree width
- Abo Khamis et al. 2016 — FAQ framework
- Wang et al. 2020 — SPORES (RA complete for LA optimization)
- Aznaveh et al. 2020, Wolf et al. 2017 — GraphBLAS triangle counting

## Git

- Remote: `git@github-rai:lums-rai/sparseblas-vs-wcoj-tc.git`
- Uses `github-rai` SSH alias

## Style Preferences

- No emojis in documents
- Show side-by-side comparisons (fused join <-> CSR code)
- Use N(x) notation, not trie operations
- "Fused" not "trie-join" when describing the WCOJ approach for HPC audience
