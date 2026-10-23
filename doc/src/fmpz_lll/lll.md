---
title: Reducing a lattice basis
functions: fmpz_lll fmpz_lll_with_removal
---

`fmpz_lll` replaces a basis of an integer lattice by an LLL-reduced
basis of the same lattice, whose first vector is within a factor
exponential in the dimension of a shortest nonzero vector and in
practice usually much closer. `fmpz_lll_with_removal` does the same and
also reports how many leading vectors have small Gram–Schmidt norms,
which is what algorithms that search for short vectors, such as van
Hoeij's factoring algorithm, need. Both choose automatically between
floating-point arithmetic of different precisions and certify the result.

::: summary algorithms/linear-algebra/lll :::

## Usage {#usage}

```c
fmpz_lll_t fl;
fmpz_lll_context_init_default(fl);  /* delta = 0.99, eta = 0.51 */
fmpz_lll(B, NULL, fl);              /* reduce the rows of B in place */
```

- *Input.* With `fl->rt == Z_BASIS`, the rows of the $d \times n$ matrix
  $B$ generate the lattice; they may be linearly dependent, in which case
  the zero vectors of the result come first, followed by a reduced basis.
  With `fl->rt == GRAM`, $B$ is the $d \times d$ Gram matrix of a basis,
  which is replaced by the Gram matrix of the reduced basis.
- *Parameters.* `fmpz_lll_context_init(fl, delta, eta, rt, gt)` sets
  $\delta \in (1/4, 1)$ and $\eta \in [1/2, \sqrt{\delta})$. Larger
  $\delta$ gives a better basis and a slower reduction. The default
  $\delta = 0.99$, $\eta = 0.51$ is close to the best possible. The Gram
  type `gt` (`APPROX` or `EXACT`) chooses whether the Gram matrix of a
  `Z_BASIS` input is maintained approximately or exactly during the
  reduction; `APPROX` is faster and the result is certified either way.
- *Transformation.* If $U$ is not `NULL`, every row operation on $B$ is
  also applied to $U$, which must have $d$ rows. Initialise $U$ to the
  identity matrix to obtain the unimodular $U$ with
  $B_{\text{new}} = U B_{\text{old}}$; a zero-initialised $U$ stays zero.
- *Output guarantee.* The result is $(\delta, \eta)$-LLL-reduced, as
  checked by `fmpz_mat_is_reduced`; in particular
  $\lVert b_1 \rVert \le (\delta - \eta^2)^{-(d-1)/2} \lambda_1$
  ([](algorithms/linear-algebra/lll#thm-lll-short)).
- *Removals.* `fmpz_lll_with_removal(B, U, gs_B, fl)` returns the number
  `newd` of leading rows to keep: the trailing rows whose squared
  Gram–Schmidt norm exceeds `gs_B` are excluded, and the first `newd` rows
  are reduced. Vectors of norm at most $\sqrt{\mathit{gs\_B}}$ in the
  lattice are combinations of the first `newd` rows. With `gs_B` equal to
  `NULL` it reduces like `fmpz_lll` and returns 0.

## Algorithm {#algorithm}

`fmpz_lll` and `fmpz_lll_with_removal` call
`fmpz_lll_with_removal_ulll` with a truncation size of 250 bits:

1. While the entries have more than 250 bits, the basis truncated to its
   leading 250 bits, with an identity matrix adjoined, is reduced; the
   transformation found in the identity block is applied to the full
   basis. This moves most of the work to small entries.
2. The full basis is reduced by `fmpz_lll_wrapper_with_removal_knapsack`:
   first L² in `double` precision (`fmpz_lll_d_with_removal_knapsack`),
   then, if that fails or the result cannot be certified, with heuristic
   exact dot products (`fmpz_lll_d_heuristic_with_removal`), and finally
   with multiprecision floating point of increasing precision
   (`fmpz_lll_mpf_with_removal`).
3. Each stage's result is certified with `fmpz_lll_is_reduced_with_removal`,
   using Villard's floating-point certificate and, if that is
   inconclusive, an exact test.

The individual stages are available as separate functions; see
[](fmpz_lll#strategies) and [](algorithms/linear-algebra/lll).

## Example {#example}

Finding the minimal polynomial of $\alpha = \sqrt{2} + \sqrt{3}$ from an
approximation: the lattice spanned by the rows $(e_i, \lfloor 2^{200}
\alpha^i \rceil)$, $i = 0, \ldots, 4$, contains the short vector
$(c_0, \ldots, c_4, \varepsilon)$ for the integer relation
$\sum c_i \alpha^i = 0$, and LLL finds it:

```c example name=fmpz_lll_minpoly
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_lll.h>

int main(void)
{
    fmpz_mat_t B, U, T;
    fmpz_lll_t fl;
    fmpz_t x, p, t, bound;
    slong i, j, d = 5, prec = 200, newd;

    fmpz_init(x);
    fmpz_init(p);
    fmpz_init(t);
    fmpz_init(bound);
    fmpz_mat_init(B, d, d + 1);
    fmpz_mat_init(U, d, d);
    fmpz_mat_init(T, d, d + 1);

    /* x = floor(2^prec sqrt(2)) + floor(2^prec sqrt(3)) ~ 2^prec alpha */
    fmpz_one(t);
    fmpz_mul_2exp(t, t, 2 * prec);
    fmpz_mul_ui(x, t, 2);
    fmpz_sqrt(x, x);
    fmpz_mul_ui(t, t, 3);
    fmpz_sqrt(t, t);
    fmpz_add(x, x, t);

    /* row i: unit vector e_i, then p ~ 2^prec alpha^i */
    fmpz_one(p);
    fmpz_mul_2exp(p, p, prec);
    for (i = 0; i < d; i++)
    {
        fmpz_one(fmpz_mat_entry(B, i, i));
        fmpz_set(fmpz_mat_entry(B, i, d), p);
        fmpz_mul(p, p, x);
        fmpz_fdiv_q_2exp(p, p, prec);
    }
    fmpz_mat_set(T, B);

    fmpz_lll_context_init_default(fl);
    fmpz_mat_one(U);
    fmpz_lll(B, U, fl);

    flint_printf("relation:");
    for (j = 0; j < d; j++)
        flint_printf(" %{fmpz}", fmpz_mat_entry(B, 0, j));
    flint_printf("  (error %{fmpz})\n", fmpz_mat_entry(B, 0, d));

    flint_printf("reduced: %d\n", fmpz_mat_is_reduced(B, fl->delta, fl->eta));
    fmpz_mat_mul(T, U, T);
    flint_printf("U B_old = B_new: %d\n", fmpz_mat_equal(T, B));

    /* only vectors of squared norm <= 1000 are wanted */
    fmpz_set_ui(bound, 1000);
    newd = fmpz_lll_with_removal(B, NULL, bound, fl);
    flint_printf("vectors kept: %wd\n", newd);

    fmpz_mat_clear(B);
    fmpz_mat_clear(U);
    fmpz_mat_clear(T);
    fmpz_clear(x);
    fmpz_clear(p);
    fmpz_clear(t);
    fmpz_clear(bound);
    return 0;
}
```
```output
relation: 1 0 -10 0 1  (error -22)
reduced: 1
U B_old = B_new: 1
vectors kept: 1
```

The relation $1 - 10\alpha^2 + \alpha^4 = 0$ gives the minimal polynomial
$x^4 - 10x^2 + 1$; the last entry, $-22$, is the rounding error of the
approximations scaled by $2^{200}$. All other vectors of the reduced basis
are far longer, so only one vector is kept with the bound 1000.
