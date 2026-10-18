---
title: Matrix multiplication over Z/nZ for multiword n
functions: mpn_mod_mat_mul
---

`mpn_mod_mat_mul` is the matrix multiplication method of the `mpn_mod`
ring: `gr_mat_mul` over the ring calls it, and so do the generic block
algorithms built on matrix multiplication (recursive LU decomposition,
triangular solving, inversion, determinants). It chooses between three
algorithms, all of which return the exact product modulo $n$:

| Algorithm | Function | Multiplications for $m \times k$ times $k \times p$ | Reductions |
|---|---|---|---|
| Classical | `gr_mat_mul_classical` with `_mpn_mod_vec_dot` | $mkp$ products of $\ell$ limbs | one per entry |
| Waksman | `mpn_mod_mat_mul_waksman` | about $(mp + m + p) k/2$ products of $\ell$ or $\ell + 1$ limbs | one per entry |
| Multimodular | `mpn_mod_mat_mul_multi_mod` | $r$ products over `nmod_mat` with $r \approx (2b + \log_2 k)/59$ | Chinese remaindering, one reduction per entry |

Here $b$ is the bit size of $n$ and $\ell$ its number of limbs.

## Choice of algorithm {#choice}

Matrices with at most 3 rows use classical multiplication. Otherwise
two cutoffs $c_W \le c_M$ are read from a table indexed by $b$ in steps
of 16 bits, tuned on x86-64:

- if $m$, $k$ or $p$ is less than $c_W$, the classical algorithm is used;
- otherwise, if one of them is less than $c_M$, Waksman's algorithm;
- otherwise the multimodular algorithm.

For moduli of at most 128 bits and of 177 to 192 bits the two cutoffs
are equal, so Waksman's algorithm is not used: the classical algorithm
runs up to dimension 88 to 230 (for $b \le 128$) or 74 (for
$177 \le b \le 192$). For other sizes Waksman's algorithm starts at
dimension 4 to 17 and the multimodular algorithm takes over between
dimensions 39 and 93; at 256 bits the cutoffs are 17 and 63, at 1024 bits
4 and 53. Strassen's algorithm is not used: at these sizes the
multimodular algorithm wins first.

## Classical multiplication {#classical}

Each entry of $C$ is a dot product of a row of $A$ and a column of $B$,
computed by `_mpn_mod_vec_dot_strided` without intermediate reductions:
the full products are added in
$\lceil (2b + \lfloor \log_2 k \rfloor + 1)/\bits \rceil$ limbs, and the
sum is reduced once ([](algorithms/arithmetic/multiword-delayed-reduction#dot)).

## Waksman's algorithm {#waksman}

Winograd [@Win1968] observed that for two vectors of even length $k$,

$$
\sum_{j=1}^{k} x_j y_j = \sum_{i=1}^{k/2} (x_{2i-1} + y_{2i})(x_{2i} + y_{2i-1}) - \sum_{i=1}^{k/2} x_{2i-1} x_{2i} - \sum_{i=1}^{k/2} y_{2i-1} y_{2i}.
$$

In a matrix product, $x$ is a row of $A$ and $y$ a column of $B$, so the
last two sums depend only on the row or only on the column and can be
computed once for each: the product costs $mkp/2 + (m + p) k/2$
multiplications instead of $mkp$. Waksman [@Waksman1970] noticed that
the correction terms can be obtained from the same kind of products.
Writing $R_l = \sum_i a_{l,2i-1} a_{l,2i}$ and
$S_t = \sum_i b_{2i,t} b_{2i-1,t}$, the identity

$$
(\alpha + \beta)(\alpha' + \beta') + (\alpha - \beta)(\alpha' - \beta') = 2(\alpha\alpha' + \beta\beta')
$$

gives $R_0 + S_t$ for every column $t$ from the products for row 0, and
$R_l + S_0$ for every row $l$ from those for column 0; then
$R_l + S_t = (R_l + S_0) + (R_0 + S_t) - (R_0 + S_0)$. The extra products
$(\alpha - \beta)(\alpha' - \beta')$ are only needed in the first row and
column. For odd $k$ the last term is added with ordinary products.

`mpn_mod_mat_mul_waksman` evaluates all sums exactly on integers and
reduces each entry of $C$ once. Accumulators have $2\ell + 1$ limbs and
hold signed values. The sums $\alpha + \beta$ of two residues fit in
$\ell$ limbs unless the top bit of $n$ is set, in which case they get an
extra limb and the products are $(\ell + 1) \times (\ell + 1)$-limb
products; the differences $\alpha - \beta$ are formed as sign and
magnitude. The halvings are exact divisions by 2 of signed multi-limb
integers. Since a multiplication of $\ell$-limb numbers costs about
$\ell^2$ word operations and an addition about $\ell$, halving the number
of multiplications at the price of $O(mkp)$ additions pays off from small
dimensions on.

## Multimodular multiplication {#multimodular}

`mpn_mod_mat_mul_multi_mod` computes the integer product $AB$ of the
matrices of residues, whose entries are less than $k (n-1)^2 < 2^{c}$
with $c = 2b + \lfloor \log_2 k \rfloor + 1$, modulo $r$ word-size moduli
whose product exceeds $2^c$, and reconstructs it by the Chinese remainder
theorem:

1. The moduli are $2^{\beta}$ and the $r - 1$ smallest primes above
   $2^{\beta}$, where $\beta$ = `NMOD_MAT_OPTIMAL_MODULUS_BITS` (59 on
   64-bit machines), with $r = 1 + \lceil (c - \beta)/\beta \rceil$. Using
   a power of two as the first modulus lets the fixed-length code of
   `flint_mpn_multi_crt` apply for small products.
2. Each row of $A$ and of $B$ is reduced modulo all moduli at once with
   `flint_mpn_multi_mod_vec`.
3. The $r$ products are computed with `nmod_mat_mul`, which itself uses
   delayed reduction, Strassen–Winograd, SIMD or BLAS kernels
   ([](algorithms/linear-algebra/matrix-multiplication)).
4. Each row of $C$ is reconstructed with `flint_mpn_multi_crt_vec` and
   each entry reduced modulo $n$ with `mpn_mod_set_mpn`.

The reductions and the reconstruction are split over threads (by rows,
when the dimensions are in the hundreds or larger), and `nmod_mat_mul`
is multithreaded itself. The code follows the multimodular
multiplication of `fmpz_mat`
([](algorithms/linear-algebra/multimodular)).

## Example {#example}

The program below multiplies two $30 \times 30$ matrices modulo
$2^{255} - 19$ with the three algorithms and checks that they agree.

```c example name=mpn_mod_mat_mul
#include <flint/fmpz.h>
#include <flint/gr.h>
#include <flint/gr_mat.h>
#include <flint/mpn_mod.h>

int main(void)
{
    gr_ctx_t ctx;
    fmpz_t p;
    gr_mat_t A, B, C1, C2, C3;
    slong i, j, d = 30;
    int status = GR_SUCCESS;

    /* p = 2^255 - 19 */
    fmpz_init(p);
    fmpz_ui_pow_ui(p, 2, 255);
    fmpz_sub_ui(p, p, 19);
    GR_MUST_SUCCEED(gr_ctx_init_mpn_mod(ctx, p));

    gr_mat_init(A, d, d, ctx);
    gr_mat_init(B, d, d, ctx);
    gr_mat_init(C1, d, d, ctx);
    gr_mat_init(C2, d, d, ctx);
    gr_mat_init(C3, d, d, ctx);

    /* A[i][j] = ij - 7, B[i][j] = A[i][j]^40 */
    for (i = 0; i < d; i++)
        for (j = 0; j < d; j++)
        {
            status |= gr_set_si(gr_mat_entry_ptr(A, i, j, ctx), i * j - 7, ctx);
            status |= gr_pow_ui(gr_mat_entry_ptr(B, i, j, ctx),
                                gr_mat_entry_ptr(A, i, j, ctx), 40, ctx);
        }

    status |= gr_mat_mul_classical(C1, A, B, ctx);
    status |= mpn_mod_mat_mul_waksman(C2, A, B, ctx);
    status |= mpn_mod_mat_mul_multi_mod(C3, A, B, ctx);

    printf("waksman: %d, multi_mod: %d\n",
        gr_mat_equal(C1, C2, ctx) == T_TRUE, gr_mat_equal(C1, C3, ctx) == T_TRUE);
    status |= gr_println(gr_mat_entry_ptr(C1, 1, 2, ctx), ctx);

    gr_mat_clear(A, ctx);
    gr_mat_clear(B, ctx);
    gr_mat_clear(C1, ctx);
    gr_mat_clear(C2, ctx);
    gr_mat_clear(C3, ctx);
    gr_ctx_clear(ctx);
    fmpz_clear(p);
    printf("status %d\n", status);
    return 0;
}
```

```output
waksman: 1, multi_mod: 1
5452489965307334725493791908297774450834361894610887641359429602394561
status 0
```

::: summary algorithms/arithmetic/multiword-delayed-reduction :::
