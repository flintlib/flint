---
title: Determinants and permanents
functions: gr_mat_det gr_mat_permanent
---

The determinant of a square matrix over a generic ring can be computed
without any division by Berkowitz's algorithm or cofactor expansion,
or with divisions by Gaussian elimination over a field or by
fraction-free elimination over an integral domain. `gr_mat_det` calls
the ring's method, which chooses among these or uses the determinant of
a specialised module. The permanent has no polynomial-time algorithm,
and `gr_mat_permanent` uses exponential-time formulas.

::: summary algorithms/linear-algebra/charpoly :::

## The determinant algorithms {#algorithms}

| Function | Algorithm | Operations | Requirements |
|---|---|---|---|
| `gr_mat_det_cofactor` | explicit formulas | $O(1)$, for $n \le 4$ | none; `GR_UNABLE` for $n > 4$ |
| `gr_mat_det_berkowitz` | constant term of the characteristic polynomial by Berkowitz's algorithm [@Ber1984] | $O(n^4)$ | none (division-free) |
| `gr_mat_det_lu` | product of the pivots of `gr_mat_lu` | $O(n^\omega)$ | invertible pivots: a field |
| `gr_mat_det_fflu` | last pivot of `gr_mat_fflu` [@Bar1968] | $O(n^3)$ | exact division: an integral domain |

The division-free algorithms work over every commutative ring,
including rings with zero divisors such as $\Z/6\Z$, and need no zero
tests, so over inexact rings they never fail for lack of a decision.
The elimination algorithms are faster but need zero tests to find
pivots; they return `GR_UNABLE` when a pivot cannot be found or inverted
(`gr_mat_det_lu`) or a division fails (`gr_mat_det_fflu`). Over a field
the LU algorithm uses fast matrix multiplication through the recursive
decomposition. Over $\Z$ the fraction-free algorithm keeps all entries
bounded by Hadamard's bound, $|\det A| \le \prod_i \|a_i\|_2$, while
Gaussian elimination over $\mathbb{Q}$ would create large fractions.

## Algorithm selection {#selection}

`gr_mat_det` calls the ring method `GR_METHOD_MAT_DET`. The default,
`gr_mat_det_generic`, uses cofactor expansion for $n \le 4$ and
Berkowitz's algorithm otherwise, which is correct over any commutative
ring. `gr_mat_det_generic_field` (cofactors, then LU) and
`gr_mat_det_generic_integral_domain` (cofactors, then fraction-free LU)
are the corresponding choices for fields and integral domains; rings
can install them as their method. The rings with their own method are:

| Ring | Method |
|---|---|
| `fmpz`, `fmpq` | `fmpz_mat_det` (multimodular with certification, [](algorithms/linear-algebra/multimodular)), `fmpq_mat_det` |
| `nmod` | `nmod_mat_det` (elimination, see [](nmod_mat/det)), for every modulus |
| `nmod8`, `nmod32`, `fmpz_mod` | cofactors for $n \le 4$, then LU, then Berkowitz if a pivot is not invertible (composite modulus) |
| `mpn_mod` | cofactors for $n \le 4$, Berkowitz for $n = 5$, then LU with Berkowitz as fallback |
| `fq`, `fq_nmod`, `fq_zech` | `gr_mat_det_generic_field` |
| `arb`, `acb` | `arb_mat_det`, `acb_mat_det` (Gaussian elimination in ball arithmetic, with preconditioning for large matrices) |
| `ca` | `ca_mat_det` |
| `arf`, `acf`, `nfloat`, `nfloat_complex` | `gr_mat_det_generic_field` |

Other rings use the generic Berkowitz algorithm for $n \ge 5$; over a
field, `gr_mat_det_lu` is much faster for large $n$.

## Permanents {#permanent}

The permanent $\operatorname{perm}(A) = \sum_{\sigma \in S_n} \prod_i
a_{i,\sigma(i)}$ differs from the determinant only by the missing signs,
but computing it is #P-complete [@Val1979]. The algorithms are:

- `gr_mat_permanent_cofactor`: recursive expansion along a column,
  with $O(n!)$ operations, used for $n \le 4$;
- `gr_mat_permanent_ryser`: Ryser's formula [@Rys1963]
  $\operatorname{perm}(A) = (-1)^n \sum_{S \subseteq \{1, \ldots, n\}}
  (-1)^{|S|} \prod_i \sum_{j \in S} a_{ij}$, with the subsets
  traversed in Gray code order so that each step updates the row sums
  by one column, $O(2^n n)$ operations;
- `gr_mat_permanent_glynn`: Glynn's formula [@Gly2010], equivalent to a
  method of Nijenhuis and Wilf [@NW1978], which sums over the $2^{n-1}$
  sign vectors $\delta$ with $\delta_1 = 1$ the products $\prod_k \delta_k
  \prod_i \sum_j \delta_j a_{ij}$ and divides by $2^{n-1}$: half the
  operations of Ryser's formula, but it needs exact division by 2;
- `gr_mat_permanent_glynn_threaded`: Glynn's formula with the sum split
  over threads.

The default `gr_mat_permanent_generic` uses cofactors for $n \le 4$, and
otherwise Glynn's formula if division by 2 may be possible (threaded
for $n \ge 11$ when threads are available and the ring is thread-safe),
falling back to Ryser's formula if that fails. `fmpz` and `fmpq` use
`fmpz_mat_permanent` and the corresponding rational function.

## Example {#example}

The program computes the determinant of a $5 \times 5$ integer matrix
reduced modulo 6, where the ring has zero divisors: the default
(Berkowitz) and the explicitly division-free algorithm succeed and agree
with the integer determinant $187 \equiv 1 \pmod 6$, while LU
decomposition stops at the pivot 2, which is not invertible. It then
computes a determinant and a permanent over $\Z[x]$.

```c example name=gr_mat_det_rings
#include <flint/flint.h>
#include <flint/gr.h>
#include <flint/gr_mat.h>

int main(void)
{
    gr_ctx_t R, Zx;
    gr_mat_t A, P;
    gr_ptr d;
    int status;
    const char * s = "[[2, 1, 0, 0, 1], [1, 3, 1, 0, 0], [0, 1, 4, 1, 0], "
                     "[0, 0, 1, 5, 1], [1, 0, 0, 1, 3]]";

    /* the same 5 x 5 matrix over Z/6Z, which has zero divisors */
    gr_ctx_init_nmod(R, 6);
    gr_mat_init(A, 5, 5, R);
    d = gr_heap_init(R);
    GR_MUST_SUCCEED(gr_mat_set_str(A, s, 0, R));

    status = gr_mat_det(d, A, R);
    flint_printf("det:           status %d, value ", status);
    gr_print(d, R); flint_printf("\n");
    status = gr_mat_det_berkowitz(d, A, R);
    flint_printf("det_berkowitz: status %d, value ", status);
    gr_print(d, R); flint_printf("\n");
    status = gr_mat_det_lu(d, A, R);
    flint_printf("det_lu:        status %d\n", status);
    gr_heap_clear(d, R);
    gr_mat_clear(A, R);
    gr_ctx_clear(R);

    /* a characteristic polynomial as a determinant over Z[x] */
    gr_ctx_init_fmpz_poly(Zx);
    gr_mat_init(P, 2, 2, Zx);
    d = gr_heap_init(Zx);
    GR_MUST_SUCCEED(gr_mat_set_str(P, "[[x - 1, -2], [-3, x - 4]]", 0, Zx));
    GR_MUST_SUCCEED(gr_mat_det(d, P, Zx));
    flint_printf("det over Z[x]: "); gr_print(d, Zx); flint_printf("\n");
    GR_MUST_SUCCEED(gr_mat_permanent(d, P, Zx));
    flint_printf("permanent:     "); gr_print(d, Zx); flint_printf("\n");
    gr_heap_clear(d, Zx);
    gr_mat_clear(P, Zx);
    gr_ctx_clear(Zx);
    return 0;
}
```
```output
det:           status 0, value 1
det_berkowitz: status 0, value 1
det_lu:        status 3
det over Z[x]: x^2-5*x-2
permanent:     x^2-5*x+10
```

The failed LU decomposition returns `GR_DOMAIN | GR_UNABLE` (3): the
inverse of the pivot does not exist, so the determinant cannot be
computed this way although it exists.
