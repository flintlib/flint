---
title: Multiplying matrices over generic rings
functions: gr_mat_mul gr_mat_mul_classical gr_mat_mul_strassen gr_mat_mul_waksman
---

`gr_mat_mul` multiplies two matrices over any ring by calling the ring's
matrix multiplication method. The other functions below
implement one algorithm each: the classical algorithm with dot
products, one level of Strassen-Winograd, and the commutative
algorithms of Waksman and Rosowski, which trade multiplications for
additions. Since block LU decomposition, triangular solving, rank,
nullspace, inverse and characteristic polynomial algorithms all call
`gr_mat_mul`, the method chosen here determines their speed as well.

::: summary algorithms/linear-algebra/overview :::

## Algorithm selection {#selection}

The default method `gr_mat_mul_generic` is classical multiplication,
which is all that can be assumed of an arbitrary ring. The rings below
override the method `GR_METHOD_MAT_MUL`:

| Ring | Method |
|---|---|
| `nmod` | `nmod_mat_mul` ([](algorithms/linear-algebra/matrix-multiplication)) |
| `nmod8` | `nmod_mat_mul_u8` from dimension 8, classical below |
| `nmod32` | the integer SIMD kernel of `nmod_mat` from dimension 8, with one level of Strassen-Winograd on top of it (single-threaded) from an architecture-dependent dimension; without the kernel, Strassen-Winograd from dimension 256 and classical below |
| `fmpz_mod` | `fmpz_mat_mul` followed by a reduction of the entries |
| `mpn_mod` | classical, Waksman, then multimodular reduction to `nmod_mat` products, with cutoffs tabulated by modulus size |
| `fmpz`, `fmpq` | `fmpz_mat_mul` ([](algorithms/linear-algebra/integer-matrix-multiplication)), `fmpq_mat_mul` |
| `fq`, `fq_nmod`, `fq_zech` | `fq_mat_mul` and its variants (Kronecker substitution into a matrix of polynomials or integers) |
| `arb`, `acb` | `arb_mat_mul`, `acb_mat_mul` (block algorithm with exact integer matrix products and error bounds, [](algorithms/ball-arithmetic/matrix-multiplication)) |
| `arf`, `acf` | classical below dimension 40 to 120 (depending on the precision), the ball algorithm above |
| `nfloat`, `nfloat_complex` | classical, fixed-point and block algorithms tuned by precision |
| `ca` | `ca_mat_mul` |

Rings without an override, such as polynomial rings, number fields and
`qqbar`, use the classical algorithm. To use Strassen or Waksman
multiplication recursively for such a ring, its context must be given
a method that calls `gr_mat_mul_strassen` (or the other algorithm)
above a cutoff, since these functions perform one level of recursion
and compute their sub-products with `gr_mat_mul`.

## The algorithms {#algorithms}

`gr_mat_mul_classical` computes each of the $mp$ entries of an $m \times
n$ by $n \times p$ product as one dot product of length $n$ with
`_gr_vec_dot`, after transposing the right factor shallowly
(`GR_MAT_TMP_INIT_SHALLOW_TRANSPOSE`) so that its columns are
contiguous. Rings with fast dot products (delayed modular reduction,
exact accumulation of integer products, a single rounding for floating
point) therefore make the classical algorithm fast as well.

`gr_mat_mul_strassen` splits each factor into $2 \times 2$ blocks and
computes the product with 7 block products and 15 block additions in
Bodrato's ordering of Winograd's variant [@Win1971; @Bodrato2010];
odd dimensions are handled by peeling off the last row or column. One
level saves an eighth of the multiplications; a full recursion costs
$O(n^{\log_2 7})$ operations [@Str1969]. See
[](algorithms/linear-algebra/strassen-winograd).

`gr_mat_mul_waksman` [@Waksman1970] computes $\sum_k a_{ik} b_{kj}$ in
pairs, as $\sum_k (a_{i,2k} + b_{2k+1,j})(a_{i,2k+1} + b_{2k,j})$ minus
correction terms that depend only on $i$ or only on $j$, in the style of
Winograd's inner product algorithm [@Win1968]. This uses about $mnp/2$
multiplications instead of $mnp$, but more additions and a division by 2
(`gr_mul_2exp_si`); it requires a commutative ring. It pays off when
multiplications are much more expensive than additions: integers and
polynomials with large entries, multi-word modular integers (as in
`mpn_mod`). `gr_mat_mul_rosowski` [@Rosowski2023] achieves the same
saving without division by 2, and saves more when the inner dimension
$n$ is odd.

All these functions return `GR_DOMAIN` when the shapes do not match.
`gr_mat_mul` passes its arguments to the ring method without checking
them, so the check is part of the contract of every `GR_METHOD_MAT_MUL`
implementation; the ring-specific methods in FLINT make it before
calling the type's own multiplication. (In FLINT 3.6 and
earlier most of them did not, and incompatible shapes gave wrong results
or corrupted memory.)

## Example {#example}

All algorithms give the same product. The program multiplies random
integer matrices of shapes $30 \times 40$ and $40 \times 35$ with the
ring method (here `fmpz_mat_mul`) and with each generic algorithm.

```c example name=gr_mat_mul_algorithms
#include <flint/flint.h>
#include <flint/gr.h>
#include <flint/gr_mat.h>

int main(void)
{
    gr_ctx_t ZZ;
    gr_mat_t A, B, C, D;
    flint_rand_t state;

    flint_rand_init(state);
    gr_ctx_init_fmpz(ZZ);
    gr_mat_init(A, 30, 40, ZZ);
    gr_mat_init(B, 40, 35, ZZ);
    gr_mat_init(C, 30, 35, ZZ);
    gr_mat_init(D, 30, 35, ZZ);

    GR_MUST_SUCCEED(gr_mat_randtest(A, state, ZZ));
    GR_MUST_SUCCEED(gr_mat_randtest(B, state, ZZ));

    /* the ring's method: fmpz_mat_mul */
    GR_MUST_SUCCEED(gr_mat_mul(C, A, B, ZZ));

    GR_MUST_SUCCEED(gr_mat_mul_classical(D, A, B, ZZ));
    flint_printf("classical: %d\n", gr_mat_equal(C, D, ZZ) == T_TRUE);
    GR_MUST_SUCCEED(gr_mat_mul_strassen(D, A, B, ZZ));
    flint_printf("strassen:  %d\n", gr_mat_equal(C, D, ZZ) == T_TRUE);
    GR_MUST_SUCCEED(gr_mat_mul_waksman(D, A, B, ZZ));
    flint_printf("waksman:   %d\n", gr_mat_equal(C, D, ZZ) == T_TRUE);
    GR_MUST_SUCCEED(gr_mat_mul_rosowski(D, A, B, ZZ));
    flint_printf("rosowski:  %d\n", gr_mat_equal(C, D, ZZ) == T_TRUE);

    /* the generic algorithms check the shapes */
    flint_printf("A*A:       status %d\n", gr_mat_mul_classical(D, A, A, ZZ));

    gr_mat_clear(A, ZZ);
    gr_mat_clear(B, ZZ);
    gr_mat_clear(C, ZZ);
    gr_mat_clear(D, ZZ);
    gr_ctx_clear(ZZ);
    flint_rand_clear(state);
    return 0;
}
```
```output
classical: 1
strassen:  1
waksman:   1
rosowski:  1
A*A:       status 1
```
