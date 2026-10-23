---
title: Linear algebra over general and inexact rings
summary: How gr_mat functions behave over rings that are not fields and over rings with undecidable zero tests, what GR_DOMAIN and GR_UNABLE mean for ranks, solutions and determinants, and which functions to use in each case.
order: 1
---

Gaussian elimination needs two things from the ring: deciding whether a
candidate pivot is zero, and inverting (or exactly dividing by) the
pivot. The generic matrix functions work over any ring but are
explicit about these two needs, so that they return an error instead
of a wrong answer when the ring does not provide them.

## Fields with exact arithmetic {#fields}

Over $\mathbb{Q}$, finite fields, number fields and `qqbar`, every zero test is
decided and every nonzero pivot is invertible. Ranks, echelon forms,
nullspaces, inverses and solutions are then always computed, and the
status codes have their plain meaning: `GR_DOMAIN` from
`gr_mat_nonsingular_solve` or `gr_mat_inv` means that the matrix is
singular, and from `gr_mat_solve_field` that the system has no
solution. The cost is that of Gaussian elimination, $O(n^3)$ field
operations, or $O(n^\omega)$ with the recursive LU decomposition when
the ring has fast matrix multiplication. Over $\mathbb{Q}$, elimination with
fractions makes the entries grow, and the fraction-free algorithms
(`gr_mat_fflu`, `gr_mat_rank_fflu`, `gr_mat_nonsingular_solve_den`) or
the multimodular algorithms of `fmpq_mat` are preferable for large
inputs.

## Integral domains and other rings {#rings}

Over an integral domain that is not a field, such as $\Z$, $\Z[x]$ or
$\mathbb{F}_p[x]$, a pivot may be nonzero without being invertible.
The fraction-free functions divide only exactly, by Sylvester's
identity, and compute in the fraction field without leaving the ring:

- `gr_mat_det` is exact over every commutative ring, by
  division-free or fraction-free algorithms;
- `gr_mat_rank` uses fraction-free elimination, and the rank is the
  rank over the fraction field;
- `gr_mat_nonsingular_solve_den` returns $X$ and $d$ with $AX = dB$, and
  `gr_mat_nonsingular_solve` succeeds only if the solution itself has
  entries in the ring (`GR_DOMAIN` otherwise);
- `gr_mat_rref_den` and `gr_mat_nullspace` return echelon forms and
  nullspace bases scaled to have entries in the ring;
- `gr_mat_inv` returns `GR_DOMAIN` unless the matrix is invertible
  over the ring itself (determinant a unit); the inverse over the
  fraction field is given by `gr_mat_adjugate` together with the
  determinant.

Over rings with zero divisors, such as $\Z/6\Z$, rank and echelon forms
are not defined in the usual sense, and `gr_mat_rank` returns
`GR_DOMAIN` when the ring is known not to be an integral domain. The
determinant, the characteristic polynomial (Berkowitz's algorithm) and
the adjugate remain available since they need no division. The
elimination-based functions fail with `GR_DOMAIN` or `GR_UNABLE` as soon
as they meet a non-invertible pivot. A ring $\Z/n\Z$ with $n$ probably
prime can be told to pretend to be a field
(`gr_ctx_set_is_pretend_field`); the field algorithms then run, and a
failure reveals a factor of $n$ through `gr_ctx_recover_zero_divisor`.

## Inexact rings {#inexact}

Over real and complex balls (`gr_ctx_init_real_arb`,
`gr_ctx_init_complex_acb`) every entry is an enclosure of an exact
value, and every result encloses the exact result for all values in
the input enclosures. A pivot such as $[0 \pm 10^{-10}]$ cannot be
shown to be zero or nonzero. The pivot search then returns `GR_UNABLE`,
and so do the rank, the solving functions and the inverse: the matrix
might be singular. When the uncertainty does not matter, the result is
computed: a well-conditioned system is solved with error bounds, and
an exactly singular matrix (with exact zero pivots) is recognised as
singular. Division-free algorithms and the ball determinant of
`arb_mat_det` always succeed and return an enclosure, which contains
zero if the matrix may be singular.

```c example name=gr_mat_exactness_arb
#include <flint/flint.h>
#include <flint/gr.h>
#include <flint/gr_mat.h>

static void
try_matrix(const char * s, gr_ctx_t RR)
{
    gr_mat_t A, B, X;
    gr_ptr d;
    slong rank;
    int status;

    gr_mat_init(A, 2, 2, RR);
    gr_mat_init(B, 2, 1, RR);
    gr_mat_init(X, 2, 1, RR);
    d = gr_heap_init(RR);
    GR_MUST_SUCCEED(gr_mat_set_str(A, s, 0, RR));
    GR_MUST_SUCCEED(gr_mat_set_str(B, "[[1], [1]]", 0, RR));

    flint_printf("%s\n", s);
    status = gr_mat_det(d, A, RR);
    flint_printf("  det    status %d: ", status); gr_print(d, RR); flint_printf("\n");
    status = gr_mat_rank(&rank, A, RR);
    flint_printf("  rank   status %d", status);
    if (status == GR_SUCCESS)
        flint_printf(": %wd", rank);
    flint_printf("\n");
    status = gr_mat_nonsingular_solve(X, A, B, RR);
    flint_printf("  solve  status %d\n", status);

    gr_heap_clear(d, RR);
    gr_mat_clear(A, RR);
    gr_mat_clear(B, RR);
    gr_mat_clear(X, RR);
}

int main(void)
{
    gr_ctx_t RR;
    gr_ctx_init_real_arb(RR, 64);
    try_matrix("[[1, 2], [3, 4]]", RR);
    try_matrix("[[1, 2], [2, 4]]", RR);
    try_matrix("[[1, 2], [2, [4 +/- 1e-10]]]", RR);
    gr_ctx_clear(RR);
    return 0;
}
```
```output
[[1, 2], [3, 4]]
  det    status 0: -2.0000000000000000000
  rank   status 0: 2
  solve  status 0
[[1, 2], [2, 4]]
  det    status 0: 0
  rank   status 0: 1
  solve  status 1
[[1, 2], [2, [4 +/- 1e-10]]]
  det    status 0: [+/- 1.01e-10]
  rank   status 2
  solve  status 2
```

The second matrix is exactly singular: the second pivot is computed as
an exact zero, so the rank is 1 and solving fails with `GR_DOMAIN` (1).
The third matrix is singular for one value in its enclosure and
nonsingular for all others, so rank and solution are undecidable and
the functions return `GR_UNABLE` (2), while the determinant is the
enclosure $[\pm 1.01 \cdot 10^{-10}]$.

A program that wants an answer anyway has several options: increase the
precision (if the inputs are exact numbers rounded to balls), use
floating-point types (`nfloat`, `arf`), which pivot by magnitude and
treat values as exact, or use exact types such as `fmpq` or `ca`. The
general pattern of trying a fast inexact method first and falling back
to an exact one on `GR_UNABLE`, while treating `GR_DOMAIN` as a definite
answer, is described in [](algorithms/generic/overview).

The exact real and complex numbers of `ca` decide most zero tests but
not all; over `ca` the pivot search (a ring method) works harder to
decide zero tests before returning `GR_UNABLE`.
