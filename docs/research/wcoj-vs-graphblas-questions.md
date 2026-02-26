# Research Questions: When Is WCOJ Strictly Necessary?

These questions frame the broader argument of the paper beyond triangle
counting. The central thesis: there exist graph computations where the
sparse linear algebra (GraphBLAS) formulation has provably worse complexity
bounds than the standard hand-written algorithm, and relational algebra
transformations are necessary and sufficient to bridge the gap.


## The Questions

### Q1. What classes of queries require WCOJ for optimal complexity?

The database community has precisely characterized this. The dividing line
is **cyclicity of the query hypergraph**:

- **Acyclic queries** (fractional hypertree width = 1): Yannakakis's
  algorithm (1981) computes any acyclic join in O(N + |output|) time using
  semijoins. Binary join plans suffice. WCOJ provides no asymptotic benefit.

- **Cyclic queries** (fractional hypertree width > 1): Binary join-project
  plans are **provably polynomially suboptimal**. The Grohe-Marx theorem
  (SODA 2006) establishes that for any class of queries with unbounded
  fractional hypertree width, there exist instances where any binary
  join-project plan is polynomially slower than the AGM bound.

The AGM bound (Atserias-Grohe-Marx, 2013) gives the maximum output size
for a conjunctive query as N^{rho*}, where rho* is the fractional edge
cover number of the query hypergraph. WCOJ algorithms (Generic Join, NPRR
2018; Leapfrog Triejoin, Veldhuizen 2014) match this bound (up to log
factors). Binary plans achieve at best O(N^{fhtw+1}), where fhtw >=
rho*.

**Key reference**: Ngo, Re, Rudra. "Skew Strikes Back: New Developments
in the Theory of Join Algorithms." SIGMOD Record 42(4), 2013.
"Any join-project plan is destined to be slower than the best possible
run time by a polynomial factor in the data size" for cyclic queries.


### Q2. Which specific graph queries are cyclic (and thus require WCOJ)?

| Query | Graph problem | fhtw | Binary worst case | WCOJ / AGM | Gap |
|-------|--------------|------|-------------------|------------|-----|
| R(x,y) JOIN S(y,z) | 2-paths / SpMM | 1 (acyclic) | O(N + out) | O(N + out) | None |
| R(x,y) JOIN S(y,z) JOIN T(z,x) | **Triangle counting** | 3/2 | O(N^2) | O(N^{3/2}) | N^{1/2} |
| C4: 4-cycle | **4-cycle counting** | 2 | O(N^2) | O(N^{3/2}) via PANDA | N^{1/2} |
| Ck: k-cycle | **k-cycle counting** | ceil(k/2) | O(N^2) | O(N^{2-1/ceil(k/2)}) | polynomial |
| K4: 4-clique | **4-clique counting** | 2 | O(N^2) | O(N^2) | constant |
| Kk: k-clique | **k-clique counting** | k/2 | O(N^{k/2}) | O(N^{k/2}) | constant |
| Diamond (K4-e) | **Diamond counting** | 3/2 | O(N^2) | O(N^{3/2}) | N^{1/2} |
| Butterfly (K_{2,2}) | **Butterfly in bipartite** | 2 | O(N^2) | O(N^{3/2}) | N^{1/2} |
| Loomis-Whitney LW_n | (generalized) | n/(n-1) | O(N^2) for n>=3 | O(N^{n/(n-1)}) | polynomial |

The pattern: **any cyclic subgraph pattern** benefits from WCOJ. This
includes all the common graph motifs (triangles, 4-cycles, diamonds,
butterflies, k-cliques for small k, houses, wheels, etc.).

Acyclic patterns (paths, trees, stars) do NOT benefit — binary joins
are already optimal for these.


### Q3. Can these cyclic queries be expressed as GraphBLAS operations? What happens?

Yes, but the GraphBLAS formulation is forced to decompose them into
pairwise matrix operations, and this decomposition is where the complexity
blowup occurs.

**Triangle counting**: C<L> = L * L, then reduce(C). The SpGEMM L*L
enumerates all length-2 paths (wedges) — O(M^2) in the worst case — even
though only O(M^{3/2}) of them close into triangles.

**4-cycle counting**: Would require computing P = A * A (all 2-paths),
then counting entries where P[i,j] * P[j,i] > 0 (or similar). The
intermediate P can have O(N^2) entries.

**k-clique counting**: The matrix power approach A^k counts walks, not
cliques. Extracting cliques from matrix products requires multiple
pairwise operations, each potentially materializing O(N^2) intermediates.

**General pattern**: For any cyclic query expressed via matrix operations,
at least one pairwise step must produce an intermediate whose size exceeds
the AGM bound. This is not an implementation limitation — it is structural.
The pairwise algebraic framework cannot express the multi-way intersection
that WCOJ uses to avoid dead-end exploration.


### Q4. Are there graph algorithms (beyond subgraph counting) where GraphBLAS has worse complexity than standard algorithms?

Yes. The GraphBLAS community tends to present only favorable comparisons,
but several complexity gaps are documented:

#### BFS: O(m) algebraic ops vs O(n+m) standard

Standard BFS with a queue: O(n+m).

Algebraic BFS via iterated SpMV (v_{k+1} = A * v_k with masking): O(m)
total algebraic operations across all iterations. For sparse graphs where
m = O(n), this matches. But Burkhardt (2021) showed that naive algebraic
BFS can be up to 24x slower, and achieving optimal O(n) algebraic
operations requires multiplying progressively smaller submatrices — a
non-obvious optimization.

Additionally, Beamer et al. (2012) showed that direction-optimized
(push-pull) BFS dramatically improves performance on scale-free graphs.
The natural algebraic formulation (repeated SpMV) cannot express this;
Yang et al. (2018) showed how to decompose push-pull into GraphBLAS, but
it required reformulating the problem outside the natural semiring
framework.

*References*:
- Burkhardt. "Optimal Algebraic BFS for Sparse Graphs." ACM TKDD 15(5), 2021.
- Beamer, Asanovic, Patterson. "Direction-optimizing BFS." SC 2012.
- Yang et al. "Implementing Push-Pull Efficiently in GraphBLAS." ICPP 2018.

#### SSSP: Bellman-Ford O(nm) vs Dijkstra O(m + n log n)

The natural GraphBLAS formulation uses the (min,+) semiring with iterated
SpMV — this is Bellman-Ford, O(nm). Dijkstra's algorithm with a Fibonacci
heap achieves O(m + n log n), which is polynomially better for dense
graphs. Priority queues are not natural semiring primitives.

Solomonik and Besta (2023) showed a linear-algebraic Fibonacci heap is
possible but it is not standard in any GraphBLAS implementation.

*Reference*: Solomonik, Besta. "Linear-algebraic implementation of
Fibonacci heap." JGAA, 2023.

#### Connected Components: O(m log n) vs O(m alpha(n))

GraphBLAS CC algorithms (label propagation / Awerbuch-Shiloach via SpMV)
achieve O(m log n) work. Optimal Union-Find achieves O(m alpha(n)) where
alpha is the inverse Ackermann function.

*Reference*: Zhang et al. "Parallel Algorithms for Finding Connected
Components using Linear Algebra." JPDC, 2020.

#### The FOSDEM 2020 slides (Szarnyas)

Gabor Szarnyas's FOSDEM 2020 talk on GraphBLAS includes a slide (37)
comparing complexity bounds. The presentation covers cases where the
algebraic formulation introduces overhead vs specialized algorithms.

*Source*: https://archive.fosdem.org/2020/schedule/event/graphblas/


### Q5. Is relational algebra strictly necessary for the transformation from bad-bound LA to optimal hand-written code?

**Yes, and this is provable via the SPORES result.**

Wang et al. (VLDB 2020) proved that relational algebra rewrite rules are
**complete** for optimizing linear algebra expressions, while linear
algebra rewrite rules alone are **not complete** for optimizing relational
algebra expressions. Formally:

- Any linear algebra expression can be rewritten to its optimal form using
  relational algebra rewrites.
- There exist relational algebra expressions (specifically, cyclic joins)
  that cannot be optimized to their optimal form using only linear algebra
  rewrites.

The reason: relational algebra can reason through higher-arity intermediate
relations (e.g., ternary relations in the triangle query), while linear
algebra is restricted to binary (matrix) operations. The multi-way
intersection that WCOJ uses is a ternary operation on the joined
variables — it cannot be decomposed into a sequence of binary matrix
operations without the pairwise blowup.

**This directly answers the question**: to go from a GraphBLAS triangle
counting expression (which has O(M^2) worst-case bounds) to the optimal
hand-written code (which has O(M^{3/2}) bounds), you MUST pass through
relational algebra. The transformation is:

1. GraphBLAS expression: C<L> = L * L (pairwise, O(M^2))
2. Relational equivalent: R(x,y) JOIN S(y,z) JOIN T(z,x) (still pairwise if evaluated left-to-right)
3. WCOJ rewrite: variable-at-a-time evaluation with multi-way intersection (O(M^{3/2}))
4. CSR instantiation: the hand-written nested-loop + galloping-intersection code

Step 3 is a relational algebra transformation (variable reordering +
intersection introduction) that has no linear algebra analog. No sequence
of matrix identity rewrites (associativity, distributivity, masking
reordering) can produce the fused three-way intersection.

*Reference*: Wang, Hutchison, Leang, Howe, Suciu. "SPORES: Sum-Product
Optimization via Relational Equality Saturation." PVLDB 13(12), 2020.


### Q6. What graph structures trigger the worst-case blowup for pairwise approaches?

The blowup is driven by **high-degree hubs** in graphs with skewed degree
distributions:

- **Star graphs**: A hub connected to sqrt(M) leaves creates M length-2
  paths through the hub, most of which don't close into triangles.

- **Power-law graphs** (social networks, web graphs): Degree exponent
  tau < 3 causes the sum of squared degrees to grow super-linearly. The
  wedge-to-triangle ratio can be 100:1 or worse.

- **RMAT/Kronecker graphs**: Standard Graph500 benchmarks. Increasing the
  RMAT skew parameter (a) concentrates edges on fewer hubs, increasing
  the intermediate blowup.

- **Real-world examples**: com-Orkut (3M nodes, 117M edges, 627M
  triangles), com-Friendster (65M nodes, 1.8B edges, 4.2B triangles),
  Twitter follower graph. All have extreme degree skew.

The WCOJ approach is immune to this blowup because it intersects eagerly
at each level — dead-end paths are never explored.


## The Proof: Why You Must Pass Through Relational Algebra

### The question

Starting from the linear algebra expression C<L> = L * L (GraphBLAS
triangle counting), can you derive the optimal hand-written galloping
intersection algorithm using only linear algebra rewrites? Or must you
"pass through" relational algebra?

**Answer: you must pass through relational algebra.** Here is the argument
at three levels of formality, aimed at linear algebra practitioners.


### Level 1: The contraction argument (most intuitive)

Matrix multiplication **contracts** a shared index by summing over it:

    C[i,k] = Sigma_j  L[i,j] * L[j,k]

The sum over j runs over ALL j in N(i), for every output entry (i,k).
This is where the O(M^2) cost comes from: the contraction explores all
length-2 paths (i -> j -> k), including those where (i,k) is not an edge.

The hand-written code does something fundamentally different. It **binds**
j in an outer loop rather than contracting it:

    for each edge (i,j):
        for each k in N(i) INTERSECT N(j):     # k constrained by BOTH i and j
            count += 1

With j bound (not summed away), the search for k is constrained by both
N(i) and N(j) simultaneously. Dead-end paths are never explored.

**The transformation from Sigma_j (contraction) to "for each j" (binding)
is not a matrix identity.** No sequence of matrix rewrites — associativity,
distributivity, transpose, masking reorder — can turn a summation into
a loop binding. This is the core of the argument.

Contraction destroys information: after computing C = L*L, the matrix C
records HOW MANY 2-paths connect each pair, but not WHICH intermediate
vertices j were used. The mask L then filters C, but the information about
j is already gone. The optimal algorithm needs j to be alive when searching
for k — and that requires relational algebra's variable binding, not
linear algebra's index contraction.


### Level 2: The dimensionality argument (for tensor people)

Matrices are rank-2 tensors: two indices. Matrix multiplication is a
binary operation that takes two rank-2 inputs and produces a rank-2 output,
contracting one shared index.

The triangle query is inherently **rank-3**: it asks for triples (i,j,k)
satisfying three simultaneous binary constraints:

    L[i,j] = 1  AND  L[j,k] = 1  AND  L[i,k] = 1

To evaluate this using only rank-2 operations, you must project the 3D
problem into 2D at some point. Any such projection contracts one variable,
and the result can be as large as O(M^2) — even though the 3D answer set
is at most O(M^{3/2}) by the AGM bound.

Concretely, whatever pairwise decomposition you choose:
- Join R(i,j) with S(j,k) first → contracts nothing, produces all 2-paths → O(M^2)
- Join R(i,j) with T(i,k) first → produces all (i,j,k) with two edges at i → O(M^2)
- Join S(j,k) with T(i,k) first → same story → O(M^2)

Every binary-first strategy produces an intermediate of size O(M^2) in the
worst case. This is not a failure of implementation — it is structural.
Any factorization of a 3D computation through 2D intermediates incurs this
blowup when the 3D query is cyclic.

The WCOJ approach never projects to 2D. It works variable-at-a-time in the
full 3D space:
- Fix i (1D constraint)
- Fix j constrained by N(i) (2D constraint)
- Find k in N(i) INTERSECT N(j) (3D constraint, evaluated as intersection)

The intersection is a ternary operation on (i,j,k) that cannot be factored
through any matrix (binary) intermediate without blowup.

**In tensor network language**: the triangle query is a cyclic tensor
contraction. Cyclic contractions cannot be evaluated by pairwise
contractions without intermediate blowup — this is a well-known result
in tensor network theory.


### Level 3: The SPORES formal result

Wang et al. (VLDB 2020) proved this formally using equality saturation:

**Theorem 2.4 (SPORES):** Two linear algebra expressions are semantically
equivalent if and only if their relational forms can be rewritten to each
other using relational equality rules (REQ).

In other words: RA rewrite rules are **complete** for the space of
equivalent LA expressions. If there exists a better way to compute an
LA expression, RA rewrites can find it.

The converse does not hold: LA rewrite rules (associativity, distributivity,
transpose identities, etc.) are NOT complete. There exist equivalent LA
expressions that cannot be reached from each other using only LA rewrites.

The proof goes through a normal form: every LA expression translates to
a relational "sum of R-monomials," and two expressions are equivalent iff
their normal forms match. The RA rewrite rules can navigate between all
equivalent normal forms; LA rewrites cannot, because they are restricted
to rank-2 intermediates.

**What this means concretely for triangle counting:**

Starting from: C<L> = L * L; count = reduce(C)

LA rewrites can produce:
- (L * L) .* L (mask after multiply)
- L .* (L * L) (same thing, Hadamard commutes)
- Block decompositions: (L1 + L2) * L = L1*L + L2*L
- Various transpose/symmetry identities

LA rewrites CANNOT produce:
- "For each edge (i,j), intersect N(i) with N(j)"

Because the intersection formulation requires:
1. Recognizing the expression as a three-way join (lifting to RA)
2. Choosing a variable ordering (RA operation: no LA analog)
3. Identifying participating atoms at each depth (RA operation)
4. Introducing multi-way intersection (RA operation: ternary constraint)
5. Lowering back to CSR operations

Steps 2-4 are relational algebra transformations with no linear algebra
counterparts. The variable ordering is a relational concept (which variable
to bind in the outer loop); the multi-way intersection is a relational
concept (simultaneously constraining a variable from multiple atoms).


### The masked SpGEMM subtlety

An important nuance: SuiteSparse:GraphBLAS implements "masked SpGEMM"
using a dot-product method that, for each edge (i,j) in the mask, computes
the inner product of row i and column j. This effectively performs the
same neighbor-list intersection that WCOJ derives.

But this is an **implementation optimization**, not an algebraic derivation.
Davis engineered the masked dot-product SpGEMM by recognizing that the
mask allows skipping irrelevant inner products — a correctness insight,
not a rewrite rule. The GraphBLAS specification says "compute C = A*B
with mask M"; the dot-product method is one implementation strategy.

The point: the LA **framework** (the algebraic specification and its
rewrite rules) cannot derive the optimization. Individual implementations
can and do implement it, but they do so by reasoning outside the algebraic
framework — effectively performing the same ternary reasoning that
relational algebra formalizes.

This is exactly the SPORES claim: you can always find the optimal
implementation by going through RA. You cannot always find it by staying
within LA.


### Why masked SpGEMM doesn't generalize beyond triangles

The masked dot-product SpGEMM works for triangles because of three
coincidences specific to the triangle query:

1. **The mask is the input**: The closing relation R(x,z) is the same
   edge set L being multiplied. Known in advance, sparse (M entries).

2. **The output is edge-indexed**: Only M inner products needed (one
   per edge in the mask).

3. **One mask application suffices**: A single inner product per edge
   fuses the three-way check.

For other cyclic patterns, these break:

**4-cycles**: P = A*A gives the 2-path counts. The 4-cycle count comes
from C(P[a,c], 2) for all pairs with P[a,c] >= 2. There is NO natural
mask — you need P[a,c] for all pairs, and the set of relevant pairs
is unknown until P is computed. For star graphs, P is dense (O(n^2)),
so the "mask" doesn't help.

**Bow-ties**: Two triangles sharing a vertex. Requires computing
per-vertex triangle counts then combining — multiple stages with
compounding intermediates, no single masked SpGEMM.

**k-cycles (k > 3)**: Same as 4-cycles — intermediate path matrices
are the bottleneck, no edge-based mask applies.

**General rule**: Masked SpGEMM works when the closing constraint is
the SAME sparse edge set used in the multiplication. This holds for
triangles (closing edge ∈ L) but not for patterns where the closing
condition involves counts, pairs of paths, or constraints between
non-adjacent vertices.

This is now Section VIII of the paper (sec-beyond.tex).


## Key References (Complete List)

### WCOJ Theory
- Atserias, Grohe, Marx. "Size bounds and query plans for relational joins." SIAM J. Computing, 2013.
- Ngo, Porat, Re, Rudra. "Worst-case optimal join algorithms." JACM 65(3), 2018.
- Ngo, Re, Rudra. "Skew Strikes Back." SIGMOD Record 42(4), 2013.
- Veldhuizen. "Leapfrog Triejoin." ICDT, 2014.
- Grohe, Marx. "Constraint Solving via Fractional Edge Covers." SODA 2006 / ACM Trans. Algorithms 2014.
- Abo Khamis, Ngo, Rudra. "FAQ: Questions Asked Frequently." PODS, 2016.
- Marx. "Tractable hypergraph properties for constraint satisfaction and conjunctive queries." JACM 60(6), 2013.
- Abo Khamis, Ngo, Suciu. "PANDA: Query Evaluation in Submodular Width." arXiv:2402.02001, 2024.

### Relational vs Linear Algebra
- Wang, Hutchison, Leang, Howe, Suciu. "SPORES." PVLDB 13(12), 2020.
- Wang, Willsey, Suciu. "Free Join: Unifying Worst-Case Optimal and Traditional Joins." SIGMOD, 2023.

### GraphBLAS and Sparse LA
- Kepner et al. "Mathematical Foundations of the GraphBLAS." HPEC, 2016.
- Davis. "Algorithm 1000: SuiteSparse:GraphBLAS." ACM TOMS 45(4), 2019.
- Wolf et al. "Fast Linear Algebra-Based Triangle Counting with KokkosKernels." HPEC, 2017.
- Aznaveh et al. "Parallel Triangle Counting and Enumeration Using Matrix Algebra." IPDPSW, 2020.
- Burkhardt. "Optimal Algebraic BFS for Sparse Graphs." ACM TKDD 15(5), 2021.
- Yang et al. "Implementing Push-Pull Efficiently in GraphBLAS." ICPP, 2018.
- Solomonik, Besta. "Linear-algebraic implementation of Fibonacci heap." JGAA, 2023.
- Zhang et al. "Parallel Algorithms for Finding Connected Components using Linear Algebra." JPDC, 2020.
- Szarnyas. "GraphBLAS: A linear algebraic approach for high-performance graph algorithms." FOSDEM, 2020.

### Graph Processing Systems Using WCOJ
- Aberger, Lamb, Tu, Notlzi, Olukotun, Re. "EmptyHeaded." ACM TODS, 2017.
- Mhedhbi, Salihoglu. "Optimizing subgraph queries by combining binary and worst-case optimal joins." PVLDB 12(11), 2019.

### Triangle Counting and Graph Algorithms
- Chiba, Nishizeki. "Arboricity and Subgraph Listing Algorithms." SIAM J. Computing 14(1), 1985.
- Latapy. "Main-memory triangle computations for very large (sparse (power-law)) graphs." TCS, 2008.
- Beamer, Asanovic, Patterson. "Direction-optimizing BFS." SC, 2012.
- Tsourakakis, Kang, Miller, Faloutsos. "DOULION: Counting triangles in massive graphs with a coin." KDD, 2009.


## How This Feeds Into the Paper

The paper currently makes the argument for triangle counting specifically.
This research suggests we can strengthen the argument:

1. **Triangle counting is not a special case** — it is the simplest
   instance of a general pattern (cyclic queries) where WCOJ strictly
   dominates pairwise approaches.

2. **The paper's SPORES argument generalizes**: the incompleteness of
   LA rewrites applies to ALL cyclic queries, not just triangles.

3. **GraphBLAS has complexity gaps beyond subgraph counting**: BFS, SSSP,
   CC all have documented cases where the natural algebraic formulation
   is suboptimal. This suggests the limitation is structural, not
   query-specific.

4. **Possible extension**: a "beyond triangles" paragraph in the
   conclusion or a brief subsection showing the 4-cycle and k-clique
   cases, with the complexity table above.
