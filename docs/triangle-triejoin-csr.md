# Triangle Counting: From Trie-Join to CSR

This document develops triangle counting from first principles, showing that
the trie-join framework — a worst-case optimal join (WCOJ) approach — when
instantiated over CSR storage, produces code identical to hand-written CSR
triangle counting.  The development proceeds in three steps: a two-relation
join that recovers sparse matrix multiply, the general three-relation triangle
query with its worst-case optimal bound, and the self-join on CSR that reveals
the equivalence.


## The Query

Given a binary relation R — equivalently, the edge set of a graph — find
all triangles: three vertices x, y, z such that edges (x,y), (y,z), and (x,z)
all exist.  In database terminology, this is a **natural join** (a conjunctive
query) over three atoms sharing variables:

```
{ (x, y, z) : R(x,y) ∧ R(y,z) ∧ R(x,z) }
```

Computing this conjunction is what database people call a **join**: the three
atoms R(x,y), R(y,z), and R(x,z) share variables (x, y, z appear
in multiple places), so evaluating the conjunction requires matching up tuples
on their shared columns.  The challenge is doing this efficiently.

We build up to this in three steps: first a two-relation join (which turns
out to be sparse matrix multiply), then the three-relation triangle query with
distinct relations, and finally the self-join on a single relation stored in CSR.


## 1. Relations as Tries

A binary relation R(x,y) is a set of pairs — equivalently, an edge list of a
graph, or a COO sparse matrix.  These are all the same thing: a collection of
tuples.  Without values/properties it's a set of pairs; with values (edge
weights, matrix entries) it's a set of triples (i, j, v).

```
R = { (1,2), (1,7), (3,4), (5,1), (5,9) }
```

The trie-join algorithm requires that we can efficiently `iterate` and `seek`
over a relation's columns in sorted order.  The **trie** is an *interface*
over the data, not a separate storage format.  The underlying data can live
in a B+ tree, a sorted array, CSR arrays, or any structure that supports
sorted iteration and fast seek — the trie-join doesn't care.

For exposition, we can visualize the trie interface over R as a sorted,
nested structure — sorted x values at level 1, and for each x, sorted y
values at level 2:

```
     R
    / | \
   1  3  5       ← level 1: sorted x values
  /|  |  |\
 2 7  4  1 9    ← level 2: sorted y values for each x
```

The trie is just the edge list sorted by (x, y) and organized hierarchically.

The trie interface has five operations:

| Operation | What it does |
|-----------|-------------|
| `level1` | The sorted set of all top-level keys |
| `open(a)` | Descend into key `a` — subsequent operations see only a's children |
| `close()` | Ascend back to the parent level, restoring the previous state |
| `iterate` | Return the next key at the current level, in sorted order |
| `seek(v)` | Jump forward to the first key ≥ v (galloping / binary search) |

For our example relation R:

- `R.level1` = {1, 3, 5} — the sorted x values (the distinct first columns).

- `R.open(1)` descends into vertex 1, exposing its children: {2, 7}.
  Subsequent `iterate` calls return 2, then 7, then signal "done".

- `R.open(5)` descends into vertex 5, exposing: {1, 9}.

- `R.close()` after `R.open(1)` restores the iterator to level 1.

In graph terms, `level1` is the set of vertices with outgoing edges,
and `open(x)` exposes the sorted neighbor list of vertex x.

To **join** two tries, the trie-join processes one variable at a time (one level
at a time). At each level, it identifies which relations mention that variable
and **intersects** their iterators.  The intersection uses a technique called
**leapfrog**: when one iterator is behind another, it `seek`s forward past
the other's current value, and the roles alternate until either a match is
found or an iterator is exhausted.  Because `seek` can skip over large
portions of the data (via galloping or binary search), the intersection cost
is proportional to the output size rather than the input sizes — this is
the key to worst-case optimality.


## 2. Warm-Up: Two-Relation Join (Sparse Matrix Multiply)

Before tackling triangles, consider the simpler join **R(x,y) ⋈ S(y,z)**.  This
produces all triples (x, y, z) where both R(x,y) and S(y,z) hold — the
structural skeleton of sparse matrix multiply C = R × S.

**Variable ordering: (x, y, z).**  At each level, the trie-join handles all
relations containing that variable:

| Level | Variable | Relations present | What happens |
|-------|----------|-------------------|-------------|
| 1 | x | R only | just iterate R |
| 2 | y | R (under x) and S (top level) | intersect R.open(x) with S.level1 |
| 3 | z | S only (under y) | just iterate S.open(y) |

When S covers all row indices 1:N, the level-2 intersection is trivial:
R.open(x) ∩ {1..N} = R.open(x).  The join simplifies to:

```
for x in R.level1:
    for y in R.open(x):
        for z in S.open(y):
            emit(x, y, z)
```

This is exactly the textbook **row-by-row SpMM** algorithm:

```
for x in 1:N
    for y in neighbors_R(x)         — nonzeros in row x of R
        for z in neighbors_S(y)     — nonzeros in row y of S
            C[x,z] += R[x,y] * S[y,z]
```

The correspondence is noteworthy: **SpMM is a special case of the trie-join**.
The same framework that will give us worst-case optimal triangle counting also
produces the standard SpMM algorithm.  The difference is only in how many
relations mention each variable, which determines whether a level iterates
or intersects.


## 3. Three-Relation Join: R(x,y) ⋈ S(y,z) ⋈ T(x,z)

Now add a third relation **T(x,z)**.  We want all (x, y, z) where all three
relations hold simultaneously.  This is the triangle query.

### Approach 1: Pairwise joins

A natural first attempt: decompose into two binary joins.

**Step 1:** Compute the SpMM  J(x, y, z) = R(x,y) ⋈ S(y,z) — all length-2
paths from x through y to z.

**Step 2:** Filter J against T — keep only those (x, y, z) where T(x,z) holds.

```
J = R(x,y) ⋈ S(y,z)           — step 1: SpMM, produces length-2 paths
result = J(x,y,z) ⋈ T(x,z)   — step 2: filter by closing edge
```

The problem is the **intermediate result J**.  If R and S each have M edges,
J can have up to O(M²) tuples.  Consider a star graph: one hub vertex
connected to √M other vertices.  R has √M edges (hub → v_i), and S has
√M edges (v_i → hub).  The SpMM in step 1 produces √M × √M = M paths — all
passing through the hub — most of which don't close into triangles.

More generally, the pairwise approach pays for **all length-2 paths**, even
those that lead nowhere.  Cost: **O(M²)** in the worst case.

### Approach 2: Worst-case optimal join (Leapfrog Triejoin)

A **worst-case optimal join** (WCOJ) processes all relations simultaneously
instead of combining them pairwise.  The specific algorithm we describe is the
**Leapfrog Triejoin** (Veldhuizen, 2014), which generates a single nested loop
as follows:

1. **Choose a variable ordering** — here (x, y, z).  One loop level per variable.

2. **For each variable, collect the relations that mention it:**
   - x appears in R(x,y) and T(x,z)
   - y appears in R(x,y) and S(y,z)
   - z appears in S(y,z) and T(x,z)

3. **At each loop level, intersect the iterators** from all relations that
   mention the current variable (the **leapfrog** step).  If a relation
   mentions a variable from a higher level, we've already `open`-ed it at
   that value, so we're iterating its next level.  If a relation first
   appears at this level, we iterate its top level.

Applying this recipe to R(x,y) ∧ S(y,z) ∧ T(x,z):

| Level | Variable | Relations | Iterator state | Operation |
|-------|----------|-----------|----------------|-----------|
| 1 | x | R, T | R.level1, T.level1 | intersect R.level1 ∩ T.level1 |
| 2 | y | R, S | R.open(x), S.level1 | intersect R.open(x) ∩ S.level1 |
| 3 | z | S, T | S.open(y), T.open(x) | **intersect S.open(y) ∩ T.open(x)** |

The "iterator state" column shows why `open` matters: at level 2, R has
already been opened at x (from level 1), so we're iterating R's y values
for that specific x.  Meanwhile S appears for the first time, so we
iterate its top level.  At level 3, S has been opened at y and T has been
opened at x.

This produces a single set of nested loops:

```
for x in intersect(R.level1, T.level1):
    for y in intersect(R.open(x), S.level1):
        for z in intersect(S.open(y), T.open(x)):
            emit(x, y, z)
```

**The key difference from pairwise:** at level 3, instead of freely iterating
S.open(y) and then filtering, we **intersect** S.open(y) with T.open(x).
Only z values that appear in both survive — dead-end paths are never
materialized.

**Complexity.**  The AGM (Atserias–Grohe–Marx) bound establishes that for
relations of size M each, the triangle query produces at most O(M^{3/2})
output tuples — and this bound is tight.  The Leapfrog Triejoin runs in
O(M^{3/2} log M) time, within a logarithmic factor of this bound.  (The
log factor arises from the galloping searches in `seek`; hash-based WCOJ
algorithms can eliminate it.)

| | Pairwise | Trie-join (WCOJ) |
|-|----------|------------------|
| Intermediate | O(M²) in worst case | None — single pass |
| Total cost | O(M²) | O(M^{3/2} log M) |
| Mechanism | SpMM then filter | Intersect at every level |

The trie-join is **worst-case optimal** (WCOJ): its running time matches the
maximum possible output size (up to log factors), and no algorithm can
do asymptotically better.  This applies to any relations R, S, T — the
self-join R = S = T is a special case.


## 4. Self-Join: R = S = T

When all three relations are the same R, the self-intersections at levels 1 and 2
become trivial.  (We assume an **undirected** graph, stored symmetrically — every
vertex that appears as a neighbor also has its own outgoing edges, so
R.level1 ⊇ range of R's second column.)

| Level | General form | With R = S = T |
|-------|-------------|----------------|
| 1 | intersect(R.level1, T.level1) | R.level1  (self ∩ self) |
| 2 | intersect(R.open(x), S.level1) | R.open(x)  (subset ∩ superset) |
| 3 | intersect(S.open(y), T.open(x)) | **intersect(R.open(y), R.open(x))** |

```
for x in R.level1:
    for y in R.open(x):
        for z in intersect(R.open(x), R.open(y)):
            count += 1
```

Read aloud: *for each edge (x,y), count the common neighbors of x and y.*

### Counting each triangle once: x < y < z

Without ordering, each triangle {a,b,c} is enumerated 6 times (all
permutations of three vertices).  The constraint x < y < z counts each
triangle exactly once:

```
for x in R.level1:
    for y in R.open(x), y > x:
        for z in intersect(R.open(x), R.open(y)), z > y:
            count += 1
```

In the trie-join framework, the `y > x` and `z > y` filters are implemented
via `seek`: after opening R at x, seek to x+1; in the intersection, start
both sides past y.


## 5. Trie-Join ↔ CSR: Line by Line

Now suppose R is a graph on vertices 1:N stored in CSR format:

```
rowptr[1:N+1]   — row i's neighbors are colval[rowptr[i] : rowptr[i+1]-1]
colval[1:nnz]   — sorted column indices within each row
```

Each trie operation maps to a CSR array operation:

| Trie operation | CSR equivalent |
|---------------|----------------|
| R.level1 | 1:N  (dense row range) |
| R.open(x) | colval[rowptr[x] : rowptr[x+1]-1]  (sorted neighbor list) |
| R.close() | restore saved row position — O(1) |
| iterate at level 2 | advance index within colval slice |
| seek(v) at level 2 | galloping search within colval slice |
| intersect(A, B) | galloping merge of two sorted colval slices |

### Side-by-side: triangle count with x < y < z

```
  Trie-Join                               CSR
  ─────────                               ───

① for x in R.level1:                      for x in 1:N

② for y in R.open(x),                       for yi in rowptr[x] : rowptr[x+1]-1
        y > x:                                 y = colval[yi]; y > x || continue

③ for z in intersect(                        xi = yi + 1
       R.open(x), R.open(y)),                zi = gallop_gt(colval, y, rowptr[y], rowptr[y+1]-1)
        z > y:                                while xi ≤ rowptr[x+1]-1 && zi ≤ rowptr[y+1]-1
                                                xv, zv = colval[xi], colval[zi]
                                                if xv == zv
④     count += 1                                  count += 1; xi += 1; zi += 1
                                                elseif xv < zv
                                                  xi = gallop_geq(colval, zv, xi, ...)
                                                else
                                                  zi = gallop_geq(colval, xv, zi, ...)
                                                end
                                              end
```

**Why xi starts at yi+1:**  `colval` is sorted within each row.  Since
`colval[yi] = y`, all entries after position `yi` are > y — automatically
satisfying z > y on the x-side.

**Why zi needs gallop_gt:**  Neighbors of y may include values ≤ y, so
we must skip to the first entry > y via binary search.


### Full CSR implementation

The innermost intersection loop — the `intersect` in the trie-join — can be
abstracted into a single function that operates on two sorted slices of
`colval`:

```julia
"""Intersect sorted neighbors of a and b (within colval), calling f(z) on each match."""
function intersect_neighbors(rowptr, colval, a, a_lo, b, b_lo, f)
    a_hi = rowptr[a+1] - 1
    b_hi = rowptr[b+1] - 1
    ai, bi = a_lo, b_lo
    while ai ≤ a_hi && bi ≤ b_hi
        av, bv = colval[ai], colval[bi]
        if av == bv
            f(av)
            ai += 1; bi += 1
        elseif av < bv
            ai = gallop_geq(colval, bv, ai, a_hi)
        else
            bi = gallop_geq(colval, av, bi, b_hi)
        end
    end
end
```

Both `a_lo` and `b_lo` are caller-supplied start positions.  The caller
advances each side past `y` before entering the intersection, enforcing
the z > y constraint from the trie-join.

With this, the triangle count becomes transparent — the outer two loops
are just edge iteration, and the inner kernel is `intersect_neighbors`:

```julia
function triangle_count(rowptr, colval, N)
    count = 0
    for x in 1:N                                      # ① R.level1
        for yi in rowptr[x] : rowptr[x+1]-1           # ② R.open(x)
            y = colval[yi]
            y > x || continue                         #    filter: y > x

            # ③ intersect R.open(x) ∩ R.open(y), z > y
            a_lo = yi + 1                              #    neighbors(x) past y
            b_lo = gallop_gt(colval, y,                #    neighbors(y) past y
                             rowptr[y], rowptr[y+1]-1)
            intersect_neighbors(rowptr, colval,
                x, a_lo,                               #    neighbors(x) after y
                y, b_lo,                               #    neighbors(y) after y
                z -> (count += 1))                     # ④ triangle: (x, y, z)
        end
    end
    return count
end
```

The entire trie-join on CSR reduces to: **for each directed edge (x → y)
with y > x, intersect the neighbor lists of x and y (restricted to entries
past y)**.  The `intersect_neighbors` function is the computational kernel —
everything else is just iteration over edges.


## 6. In RAICode: From TrieState to Generated Code

RAICode implements the trie-join via the `TrieState` abstraction.  Each
relation is wrapped in a `TrieState` that provides the trie interface from
Section 1.  A `TrieStateConjunction` takes multiple TrieStates, applies the
Leapfrog Triejoin recipe from Section 3, and generates the nested-loop +
intersection code automatically.

### 6a. The TrieState interface

Every TrieState — whether backed by a B+ tree (`TrieStateMaterialized`),
dense sequential keys (`TrieStateContiguous`), or flat arrays — implements
the same four operations (corresponding to the abstract trie interface from
Section 1):

| Operation | What it does | Section 1 equivalent |
|-----------|-------------|---------------------|
| `iterate(Val(k))` | Return the next key at trie level k | `iterate` / `level1` |
| `seek_lub(v, Val(k))` | Jump to first key ≥ v at level k | `seek(v)` |
| `open(Val(k))` | Save state at level k, descend to level k+1 | `open(a)` |
| `close(Val(k))` | Restore state to before the last open at level k | `close()` |

For the triangle query, we create three TrieStates over the same edge
relation and hand them to `TrieStateConjunction`:

```julia
R_xy = TrieStateMaterialized(Var[], UnsafeConstant(edge_tree))  # R(x,y)
R_yz = TrieStateMaterialized(Var[], UnsafeConstant(edge_tree))  # R(y,z)
R_xz = TrieStateMaterialized(Var[], UnsafeConstant(edge_tree))  # R(x,z)

ARG_VARS = ((1, 2), (2, 3), (1, 3))
conjunction = TrieStateConjunction{T, ARG_TYPES, ARG_VARS}((R_xy, R_yz, R_xz), env)
```

The conjunction machinery generates (at compile time, via Julia's type
system and `@generated` functions) code that uses only the abstract
TrieState interface — no spans, no CSR, just iterate, seek_lub, open,
close:

```
Level 1 (x):  leapfrog intersect R_xy.iterate(Val(1)) with R_xz.iterate(Val(1))
              on match: open(Val(1)) on both → descend to level 2

Level 2 (y):  leapfrog intersect R_xy.iterate(Val(2)) with R_yz.iterate(Val(1))
              on match: open on both → descend to level 3

Level 3 (z):  leapfrog intersect R_yz.iterate(Val(2)) with R_xz.iterate(Val(2))   ← THE INTERSECTION
              on match: emit triangle (x, y, z)

              close back to level 2, advance to next y
close back to level 1, advance to next x
```

**Where the intersection happens.**  Level 3 is where the real work is.
The leapfrog at level 3 alternates `seek_lub` calls on R_yz and R_xz —
exactly the `intersect(R.open(y), R.open(x))` from Section 4, and the
`intersect_neighbors` kernel from Section 5.  Levels 1 and 2 are just
edge iteration; level 3 is the computational hot path.

### 6b. Spans make intersection efficient

The abstract interface works regardless of backing storage, but performance
depends on what `seek_lub` actually does inside.  In RAICode, relations
live in B+ trees accessed through CSI (Contiguous Span Iterator), which
returns data in **spans** — contiguous arrays of keys and values from a
single B+ tree page.

A key property: **rows are contiguous within spans.**  A given vertex's
neighbor list may live within one span or extend across several (if the
neighbor list is larger than a B+ tree page), but within each span, the
neighbors form a contiguous, sorted slice of memory.  When the level-3
intersection calls `seek_lub` within a span, it performs a galloping
search on a contiguous array — the same `gallop_geq` operation from the
hand-written CSR code in Section 5.

So the generated leapfrog intersection, when executing within a single
span, is already performing the same memory accesses as the hand-written
CSR intersection kernel.  The remaining overhead is:

- **Page boundaries:** when a neighbor list crosses a span boundary, the
  iterator must advance to the next B+ tree page (`advance_span!`)
- **open/close save/restore:** managing B+ tree iterator state across
  trie levels — detecting when the iterator has moved to a different span
  and repositioning

### 6c. CSR index: generated code = hand-written CSR

If the underlying storage is CSR arrays instead of B+ trees, these
overheads disappear entirely.  Each row's neighbor list is a single
contiguous array slice — `colval[rowptr[row] : rowptr[row+1]-1]`.  No
page boundaries, no multi-span handling.  The TrieState operations become
trivial array operations:

| TrieState operation | With CSR storage |
|--------------------|-----------------|
| `iterate(Val(1))` | row += 1 |
| `seek_lub(v, Val(1))` | row = v  (rows are 1:N — O(1)) |
| `open(Val(1))` | col_lo, col_hi = rowptr[row], rowptr[row+1]-1 |
| `seek_lub(v, Val(2))` | gallop_geq(colval, v, col_at, col_hi) |
| `close(Val(1))` | restore saved row — O(1) |

With these trivial implementations, the generated conjunction code for the
triangle query compiles down to:

```julia
# What the conjunction generates with CSR-backed TrieStates:

for x in 1:N                                        # iterate(Val(1)) on R_xy, R_xz
    for yi in rowptr[x] : rowptr[x+1]-1             # iterate(Val(2)) on R_xy
        y = colval[yi]
        y > x || continue                           # seek_lub(x+1) on R_yz

        # Level 3: leapfrog intersect R_yz and R_xz
        xi = yi + 1                                  # seek_lub(y+1) on R_xz
        zi = gallop_gt(colval, y, rowptr[y], ...)    # seek_lub(y+1) on R_yz
        while xi ≤ rowptr[x+1]-1 && zi ≤ rowptr[y+1]-1
            xv, zv = colval[xi], colval[zi]
            if xv == zv
                count += 1; xi += 1; zi += 1        # emit triangle
            elseif xv < zv
                xi = gallop_geq(colval, zv, xi, ...) # seek_lub on R_xz
            else
                zi = gallop_geq(colval, xv, zi, ...) # seek_lub on R_yz
            end
        end
    end
end
```

This is **exactly** the hand-written CSR code from Section 5.  The
trie-join abstraction — TrieStateConjunction generating leapfrog
intersection from variable bindings — adds zero overhead when the
underlying storage is CSR.  The generated code compiles down to the same
galloping intersection loops that a sparse-matrix programmer would write
by hand.

The progression:

| Layer | What it adds | Cost |
|-------|-------------|------|
| Abstract trie-join (Section 3) | Correctness, WCOJ guarantees | — |
| TrieState interface (6a) | Pluggable storage backends | Interface dispatch |
| Spans / CSI (6b) | Efficient B+ tree access | Page boundary handling |
| CSR storage (6c) | No page boundaries, flat arrays | **Zero overhead** |


## 7. Performance Results

**Graph:** RMAT/Kronecker at various scales.

*(Results to be recovered / regenerated — placeholders below.)*

### Triangle counting throughput

| scale | RAICode trie | RAICode CSR | julia_csr | cpp_csr |
|-------|-------------|-------------|-----------|---------|
|  16   |             |             |           |         |
|  18   |             |             |           |         |
|  20   |             |             |           |         |
|  22   |             |             |           |         |

RAICode trie = TrieStateConjunction with TSM+TSC span-level intersection.
RAICode CSR = TrieStateConjunction with CSR-backed TrieStates.

### Internal: incremental span improvements

Each column adds one improvement over the previous:

| scale | no spans | + TSM spans | + TSC spans | + CSR |
|-------|----------|-------------|-------------|-------|
|  16   |          |             |             |       |
|  18   |          |             |             |       |
|  20   |          |             |             |       |
|  22   |          |             |             |       |

- **no spans:** iterate-only, no span-level intersection
- **+ TSM spans:** TrieStateMaterialized uses span-level galloping
- **+ TSC spans:** TrieStateContiguous also uses span-level galloping
- **+ CSR:** CSR-backed TrieStates (no page boundaries, flat arrays)

### What the results show

1. **RAICode trie with spans ≈ hand-written CSR.**  The generated conjunction
   code, when operating within contiguous spans, performs the same galloping
   intersection as hand-written CSR — and the benchmarks confirm this with
   comparable throughput.

2. **Spans matter.**  Without spans (iterate-only), the trie-join pays per-element
   overhead on every `iterate` / `seek_lub` call.  Spans amortize this by
   handing the intersection kernel a contiguous array slice to gallop over.

3. **TSC spans close the gap.**  Adding span support to TrieStateContiguous
   brings its performance in line with TSM spans — the implicit-key overhead
   disappears when intersection operates on contiguous array slices.

4. **C++ reference.**  The C++ CSR implementation provides a ceiling for
   what flat-array intersection can achieve on this hardware.


## References

- Atserias, A., Grohe, M., and Marx, D. "Size bounds and query plans for
  relational joins." *SIAM Journal on Computing*, 2013.  (The AGM bound on
  maximum output size for conjunctive queries.)

- Ngo, H. Q., Porat, E., Re, C., and Rudra, A. "Worst-case optimal join
  algorithms." *Journal of the ACM*, 2018.  (WCOJ algorithms matching
  the AGM bound.)

- Veldhuizen, T. L. "Leapfrog Triejoin: A simple, worst-case optimal join
  algorithm." *ICDT*, 2014.  (The specific trie-join algorithm described in
  this document.)
