---
title: bool_mat.h
header: bool_mat.h
description: Dense matrices over the Boolean semiring, with products, powers, transitive closure and graph algorithms on adjacency matrices.
category: Matrices
order: 8
---

A `bool_mat_t` is a dense matrix over the Boolean semiring
$(\{0, 1\}, \vee, \wedge)$: sums are logical *or* and products are
logical *and*. Entries are stored as `int` values, which must be 0 or 1;
the arithmetic functions use bitwise operations and give wrong results
for other values. A square Boolean matrix is the adjacency matrix of a
directed graph, with $A_{ij} = 1$ when there is an edge from vertex $i$
to vertex $j$; then $A^k$ records which vertices are joined by a walk of
length exactly $k$. The graph functions documented at the end answer
questions about walks, cycles and strongly connected components.

`arb_mat_exp` and `acb_mat_exp` use `bool_mat_nilpotency_degree` and
`bool_mat_all_pairs_longest_walk` on the pattern of nonzero entries of
the input: when the pattern is
nilpotent, the Taylor series of the exponential is a finite sum, and its
error bound can be set to zero entry by entry.

The dimensions of a matrix are fixed at initialisation and may be zero.
Unless stated otherwise, the caller must ensure that the dimensions of
the arguments are compatible; some functions throw an exception if they
are not. Output and input matrices may be the same object.

## Types and entries {#types}

::: type
bool_mat_struct
bool_mat_t

A matrix: `entries` points to the $rc$ entries stored row by row, `rows`
to an array of pointers to the start of each row, and `r` and `c` are
the numbers of rows and columns. A `bool_mat_t` is an array of one
`bool_mat_struct`, so it is passed by reference.
:::

::: macro
bool_mat_nrows(mat)
bool_mat_ncols(mat)

The number of rows and of columns of the matrix.
:::

::: function
int bool_mat_get_entry(const bool_mat_t mat, slong i, slong j)
void bool_mat_set_entry(bool_mat_t mat, slong i, slong j, int x)

Get and set the entry in row $i$ and column $j$. The value $x$ must be 0
or 1. Indices are not checked.
:::

## Memory management {#memory}

::: function
void bool_mat_init(bool_mat_t mat, slong r, slong c)

Initialises `mat` as the $r \times c$ zero matrix. No memory is allocated
if $r$ or $c$ is zero.
:::

::: function
void bool_mat_clear(bool_mat_t mat)

Frees the memory used by `mat`.
:::

::: function
void bool_mat_swap(bool_mat_t mat1, bool_mat_t mat2)

Swaps the two matrices, including their dimensions, in constant time.
:::

## Assignment and special matrices {#assignment}

::: function
void bool_mat_set(bool_mat_t dest, const bool_mat_t src)

Copies `src` to `dest`, which must have the same dimensions.
:::

::: function
void bool_mat_zero(bool_mat_t mat)
void bool_mat_one(bool_mat_t mat)

Set `mat` to the zero matrix, or to the matrix with ones on the main
diagonal and zeros elsewhere (also for non-square matrices).
:::

::: function
void bool_mat_directed_path(bool_mat_t A)
void bool_mat_directed_cycle(bool_mat_t A)

Set the square matrix $A$ of order $n$ to the adjacency matrix of the
directed path $0 \to 1 \to \cdots \to n - 1$, with $A_{ij} = 1$ exactly
when $j = i + 1$, or of the directed cycle, with $A_{ij} = 1$ exactly
when $j \equiv i + 1 \pmod n$. Throw an exception if $A$ is not square.
:::

## Random generation {#random}

::: function
void bool_mat_randtest(bool_mat_t mat, flint_rand_t state)
void bool_mat_randtest_diagonal(bool_mat_t mat, flint_rand_t state)

Set `mat` to a random matrix, or to a random matrix that is zero off the
main diagonal. A density between 0 and 1 is chosen at random, and each
entry (each diagonal entry) is 1 with that probability.
:::

::: function
void bool_mat_randtest_nilpotent(bool_mat_t mat, flint_rand_t state)

Sets the square, nonempty matrix `mat` to a random nilpotent matrix: a
random strictly lower triangular matrix with rows and columns permuted
by a common random permutation, that is, the adjacency matrix of a random
directed acyclic graph. Throws an exception if `mat` is not square or is
empty.
:::

## Input and output {#io}

::: function
void bool_mat_print(const bool_mat_t mat)
void bool_mat_fprint(FILE * file, const bool_mat_t mat)

Print the matrix to standard output or to `file`, one row per line in the
form `[1, 0, 1]`.
:::

## Comparison and properties {#properties}

::: function
int bool_mat_equal(const bool_mat_t mat1, const bool_mat_t mat2)

Returns 1 if the matrices have the same dimensions and the same entries,
and 0 otherwise.
:::

::: function
int bool_mat_any(const bool_mat_t mat)
int bool_mat_all(const bool_mat_t mat)

Return 1 if some entry, respectively every entry, of `mat` is 1, and 0
otherwise. For an empty matrix `bool_mat_any` returns 0 and
`bool_mat_all` returns 1.
:::

::: function
int bool_mat_is_empty(const bool_mat_t mat)
int bool_mat_is_square(const bool_mat_t mat)

Return 1 if `mat` has zero rows or zero columns, respectively the same
number of rows and columns, and 0 otherwise.
:::

::: function
int bool_mat_is_diagonal(const bool_mat_t A)
int bool_mat_is_lower_triangular(const bool_mat_t A)

Return 1 if $A_{ij} = 0$ for all $i \ne j$, respectively for all $i < j$,
and 0 otherwise. The matrix need not be square.
:::

::: function
int bool_mat_is_transitive(const bool_mat_t mat)

Returns 1 if $A_{ij} \wedge A_{jk} \implies A_{ik}$ for all $i, j, k$,
that is, if $A^2 \le A$ entrywise, and 0 otherwise. Requires a square
matrix; uses $O(n^3)$ operations.
:::

::: function
int bool_mat_is_nilpotent(const bool_mat_t A)

Returns 1 if some power $A^k$ with $k \ge 1$ is zero, that is, if the
graph of $A$ has no cycles (including loops $A_{ii} = 1$), and 0
otherwise. Uses a depth-first search for cycles, in $O(n^2)$ operations.
Throws an exception if $A$ is not square. For an empty matrix it returns
0.
:::

## Arithmetic {#arithmetic}

::: function
void bool_mat_transpose(bool_mat_t dest, const bool_mat_t src)

Sets `dest` to the transpose of `src`. Throws an exception if the
dimensions of `dest` are not those of the transpose.
:::

::: function
void bool_mat_complement(bool_mat_t B, const bool_mat_t A)

Sets $B$ to the entrywise negation of $A$: $B_{ij} = \neg A_{ij}$.
:::

::: function
void bool_mat_add(bool_mat_t res, const bool_mat_t mat1, const bool_mat_t mat2)

Sets `res` to the sum of the two matrices, the entrywise *or*.
:::

::: function
void bool_mat_mul_entrywise(bool_mat_t res, const bool_mat_t mat1, const bool_mat_t mat2)

Sets `res` to the entrywise product of the two matrices, the entrywise
*and*. Throws an exception if the input dimensions differ.
:::

::: function
void bool_mat_mul(bool_mat_t res, const bool_mat_t mat1, const bool_mat_t mat2)
void bool_mat_sqr(bool_mat_t B, const bool_mat_t A)

Set `res` to the matrix product, $\mathit{res}_{ij} = \bigvee_k
(\mathit{mat1}_{ik} \wedge \mathit{mat2}_{kj})$, and $B$ to $A^2$. The
product is computed by the classical algorithm with early exit, in
$O(rck)$ operations. Throws an exception if the dimensions are
incompatible.
:::

::: function
void bool_mat_pow_ui(bool_mat_t B, const bool_mat_t A, ulong exp)

Sets $B$ to $A^{\mathit{exp}}$ by binary exponentiation, with $A^0$ the
identity matrix. Throws an exception if $A$ is not square.
:::

::: function
int bool_mat_trace(const bool_mat_t mat)

Returns the trace of the square matrix in the Boolean semiring: 1 if some
diagonal entry is 1, and 0 otherwise (also for an empty matrix). Throws
an exception if the matrix is not square.
:::

## Graph functions {#graphs}

These functions regard a square matrix $A$ of order $n$ as the adjacency
matrix of a directed graph on the vertices $0, \ldots, n - 1$. They throw
an exception if $A$ is not square. See [](bool_mat/graphs) for an
example.

::: function
void bool_mat_transitive_closure(bool_mat_t B, const bool_mat_t A)

Sets $B$ to the transitive closure $A + A^2 + A^3 + \cdots$ of $A$:
$B_{ij} = 1$ if and only if there is a walk of positive length from $i$
to $j$. Uses Warshall's algorithm [@War1962], in $O(n^3)$ operations.
:::

::: function
slong bool_mat_get_strongly_connected_components(slong * p, const bool_mat_t A)

Computes the strongly connected components of the graph of $A$ and
returns their number $k$. Sets $p_i$, for each vertex $i$, to the index
in $\{0, \ldots, k - 1\}$ of the component containing $i$. The components
form a directed acyclic graph, and they are numbered in postorder with
respect to it: if there is an edge from a vertex of component $a$ to a
vertex of a different component $b$, then $a > b$. Uses Tarjan's
algorithm [@Tar1972], in $O(n^2)$ operations.
:::

::: function
slong bool_mat_nilpotency_degree(const bool_mat_t A)

Returns the nilpotency degree of $A$, the smallest $k \ge 1$ with
$A^k = 0$, or $-1$ if $A$ is not nilpotent. For an empty matrix it
returns 0. For a nilpotent matrix the degree is one more than the length
of the longest path in the graph, which is computed by a topological
sort.
:::

::: function
slong bool_mat_all_pairs_longest_walk(fmpz_mat_t B, const bool_mat_t A)

Sets $B_{ij}$ to the length of the longest walk from vertex $i$ to vertex
$j$ in the graph of $A$, where the empty walk from a vertex to itself has
length 0. $B_{ij} = -1$ if there is no walk from $i$ to $j$, and
$B_{ij} = -2$ if there are arbitrarily long walks (a walk from $i$ to $j$
through a cycle). The matrix $B$ must be $n \times n$.

Returns $-2$ if some entry of $B$ is $-2$, and otherwise the largest
entry of $B$, which is one less than `bool_mat_nilpotency_degree`; for
an empty matrix it returns $-1$.

::: details
The function works on the condensation of the graph, the acyclic graph
of its strongly connected components. An application is to matrix power
series: if $A$ has the sparsity pattern of a real or complex matrix $M$
and $B_{ij} \ne -2$, then $(M^k)_{ij} = 0$ for all $k > B_{ij}$, so the
tail $\sum_{k > N} a_k (M^k)_{ij}$ of a power series vanishes for
$N \ge B_{ij}$.
:::
:::
