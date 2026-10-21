---
title: Multiplying integer matrices
functions: fmpz_mat_mul fmpz_mat_mul_blas fmpz_mat_mul_fft_small
---

`fmpz_mat_mul` computes $C = AB$ for an $m \times k$ matrix $A$ and a
$k \times n$ matrix $B$ with integer entries. It is the function to call:
from the dimensions and the sizes of the entries it chooses among about
ten methods, from dot products accumulated in one machine word to
floating-point `gemm` modulo small primes, multimodular multiplication
over $\Z/p\Z$ and transforms of the entries with the small-prime FFT.
All methods compute the exact product, so the result never depends on the
choice. The functions with a suffix force one method; they are useful for
testing and for callers that know their inputs better than the
heuristics do.

::: summary algorithms/linear-algebra/integer-matrix-multiplication :::

## Usage {#usage}

- $C$ must be initialised as an $m \times n$ matrix. The dimensions are
  checked only by `FLINT_ASSERT`.
- `fmpz_mat_mul` allows $C$ to be $A$ or $B$; the product is then
  computed into a temporary matrix. The forced methods do not allow
  aliasing. Windows are allowed, without partial overlap.
- `fmpz_mat_mul_blas`, `fmpz_mat_mul_fft_small` and
  `fmpz_mat_mul_fft_small_trunc` return 1 if they computed the product
  and 0 if their method does not apply, so that a caller can fall back to
  `fmpz_mat_mul`.
- Most methods use several threads when `flint_set_num_threads` allows
  it: the small-entry and two-word kernels split the rows of $A$, the
  multimodular method computes the residue products with threaded
  `nmod_mat_mul` and threads the reduction and reconstruction, and
  `flint_dgemm` is threaded.
- The square $A^2$ is computed by `fmpz_mat_sqr`, which uses Bodrato's
  squaring formulas for tiny dimensions and otherwise `fmpz_mat_mul`;
  powers by `fmpz_mat_pow`.

## How the method is chosen {#dispatch}

Let $d = \min(m, k, n)$, let $a$ and $b$ be the bit sizes of the largest
entries of $A$ and $B$ in absolute value (from `fmpz_mat_max_bits`), let
$s = 1$ if some entry is negative and $s = 0$ otherwise, and let
$c = a + b + \lceil \log_2 (k + 1) \rceil$, so that every entry of $C$ is
less than $2^c$ in absolute value. On a 64-bit machine `fmpz_mat_mul`
tests the following cases in order.

| Case | Method |
|---|---|
| $k = 0$, or $A$ or $B$ zero | $C = 0$ |
| $k = 1$; $k = 2$ | `fmpz_mat_mul_classical`; one `fmpz_fmma` per entry |
| $\min(a, b) \ge 10000$ ($d \le 5$), $\ge 7000 + 0.15 d^2$ ($d > 5$) | `fmpz_mat_mul_fft_small`, if it succeeds |
| $d > 50$ and $d > L(c)$ | `fmpz_mat_mul_blas`, if it succeeds |
| $a, b \le 62$, and $m < 9$ or $m + k < 20$ | `_fmpz_mat_mul_small_1`, `_2a` or `_2b` |
| $a, b \le 62$, $c \le 62$, $d > 1000 + 300t$ | `fmpz_mat_mul_strassen` |
| $a, b \le 62$, $c > 62$, $d > 4000 + 300t$ | `fmpz_mat_mul_multi_mod` |
| $a, b \le 62$ otherwise | `_fmpz_mat_mul_small_internal` |
| $a + s, b + s \le 128$ | `_fmpz_mat_mul_double_word_internal`, or `fmpz_mat_mul_multi_mod` for large $d$ |
| $d \ge 3 \lceil \log_2 (c + 1) \rceil$ | `fmpz_mat_mul_multi_mod` |
| $d < 20$, large balanced entries | `fmpz_mat_mul_waksman` |
| $a, b \ge 500$ and $d \ge 8$ | `fmpz_mat_mul_strassen` |
| otherwise | `fmpz_mat_mul_classical` |

Here $t$ is the number of threads, and $L(c)$ is $0$ for $c \le 53$,
1000 for $c \le 64$, 1300 for $c \le 128$, 250 for $c \le 256$ and
$200 + 8 \lceil \log_2 (c + 1) \rceil$ beyond. In the two-word case the
multimodular method is chosen when $d' > 300 + ((c - 128)/8)^2 t$, with
$d' = d$ for nonnegative inputs and $d' = 2d$ otherwise. Waksman's
algorithm is chosen for $d = 2$ when $\min(a, b) \ge 5000$ and
$\max(a, b) \le 1.1 \min(a, b)$, and for $3 \le d < 20$ when
$\max(a, b) \le 1.6 \min(a, b)$ and $\min(a, b)$ is at least 3000
($d = 3$), 1000 ($4 \le d < 12$) or 500 ($d \ge 12$).

## Methods {#methods}

### Dot products {#classical}

`fmpz_mat_mul_classical` copies the columns of $B$ into contiguous arrays
(shallow copies of the `fmpz` words, without allocation) and computes
each entry with `_fmpz_vec_dot_general`. It costs $mkn$ multiplications of
entries.

When all entries are small `fmpz` (at most 62 bits on a 64-bit machine),
the kernels `_fmpz_mat_mul_small_*` read the words directly and accumulate
each dot product in one word if $c \le 62$, in two words if
$c \le 127$, and in three words otherwise, converting to `fmpz` once per
entry. The blocked version `_fmpz_mat_mul_small_internal` transposes $B$
into words, processes blocks of 16 rows of $A$ against blocks of 128
columns of the inner dimension, and distributes rows over threads. For
entries up to two words, `_fmpz_mat_mul_double_word_internal` packs the
entries into two-word integers and accumulates the products of $2 \times
2$ words in four or five words.

### Waksman's algorithm {#waksman}

Waksman's algorithm [@Waksman1970] computes the products of pairs of
sums $(a_{i,2l} + b_{2l+1,j})(a_{i,2l+1} + b_{2l,j})$ and recovers the
dot products from them with corrections that depend only on the rows of
$A$ or the columns of $B$. It needs about $mkn/2$ multiplications of
entries instead of $mkn$, at the price of more additions, which pays off
for small matrices with large entries of similar sizes, where
multiplications dominate.

### Strassen–Winograd {#strassen}

`fmpz_mat_mul_strassen` performs one level of Strassen's algorithm in
Winograd's form, with the sequence of operations given by Bodrato
[@Bodrato2010], and computes the seven half-size products with
`fmpz_mat_mul`, which may recurse. Odd rows and columns are peeled off and
handled with matrix-vector products. This saves multiplications of large
entries at the cost of 15 matrix additions per level; see
[](algorithms/linear-algebra/strassen-winograd).

### Multimodular {#multi-mod}

`fmpz_mat_mul_multi_mod` reduces $A$ and $B$ modulo several word-size
moduli, multiplies the residue matrices with `nmod_mat_mul`, and
reconstructs $C$ by the Chinese remainder theorem in the symmetric range
when $s = 1$. The moduli are chosen so that their product exceeds
$2^{c + s}$: if $c + s \le 63$, the single modulus $2^{c + s}$ (a power of
two, for which `nmod_mat_mul` is as fast as for a prime); otherwise
$2^{63}$ and $\lceil (c + s - 63)/59 \rceil$ primes above $2^{59}$. The
reduction and reconstruction use the subproduct trees of `fmpz_comb_t`
(`fmpz_mat_multi_mod_ui_precomp`, `fmpz_mat_multi_CRT_ui_precomp`). Since
`nmod_mat_mul` runs at nearly the speed of floating-point arithmetic,
this is the method of choice for large matrices; see
[](algorithms/linear-algebra/integer-matrix-multiplication#multimodular).

### Floating point {#blas}

`fmpz_mat_mul_blas` uses the double precision matrix product
`flint_dgemm`, which calls a BLAS library if FLINT was built with one
and otherwise FLINT's own kernels ([](algorithms/linear-algebra/gemm)).
If $a + b + \lceil \log_2 (k + 1) \rceil \le 53$, one `dgemm` on the
entries converted to `double` gives the exact product. Otherwise the
entries are reduced modulo primes $p < 2 + 2\sqrt{(2^{53} - 1)/k}$ (at most
1664544 when $c > 200$), stored as 32-bit integers, and lifted to
`double` in the symmetric range, so that each residue product is an exact
`dgemm`; as many primes are used as needed to exceed $2^{c + s}$, and $C$
is reconstructed by the Chinese remainder theorem. The function returns 0
if $k$ is so large that $p$ would drop below 1000, if a dimension exceeds
`INT_MAX`, or on 32-bit machines.

### Transformed entries {#fft-small}

For entries of many thousands of bits, the cost is dominated by the
$mkn$ multiplications of large integers, each of which, done separately,
transforms both operands and the result. `fmpz_mat_mul_fft_small`
represents the entries in the transformed integer ring of `gr`
(`gr_ctx_init_transformed_mpn`), built on the small-prime FFT of
[](fft_small): each of the $mk + kn$ entries is transformed once, each
entry of $C$ is accumulated as a sum of $k$ pointwise products in the
transformed domain, and only the $mn$ results are transformed back. This
saves about two thirds of the transforms. Signs are handled by the
transformed ring.

The storage for the transformed entries is limited by
`flint_fft_small_max_transformed_ring_size`. When all of them do not fit,
the function blocks the computation: one operand resident and the other
streamed by rows or columns (no entry is transformed twice), then both
operands in blocks (the streamed side is transformed again for each
block), and as a last resort blocks of the inner dimension whose partial
results are added as integers. Squaring ($A = B$ as pointers) keeps a
single pool of transforms when everything is resident.

`fmpz_mat_mul_fft_small_trunc` returns, for each entry, the limbs of its
absolute value from limb `lo_limbs` upwards, with its sign; the dropped
low limbs perturb the result by at most one unit in its last limb. This
is what a working-precision product needs, for example in ball
arithmetic. It refuses the last blocking regime, in which truncation
errors would accumulate.

`fmpz_mat_mul_fft` is an older transform-based method using the
Schönhage–Strassen FFT of [](fft) for entries of at least 2000 bits. It
is no longer used by `fmpz_mat_mul`.

## Example {#example}

The product of a $30 \times 20$ and a $20 \times 25$ matrix with random
entries of 40 and 70 bits, computed by the default method and by forced
methods, which all agree; then a product in place:

```c example name=fmpz_mat_mul_methods
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>

int main(void)
{
    fmpz_mat_t A, B, C, D;
    flint_rand_t state;
    int ok;

    flint_rand_init(state);
    fmpz_mat_init(A, 30, 20);
    fmpz_mat_init(B, 20, 25);
    fmpz_mat_init(C, 30, 25);
    fmpz_mat_init(D, 30, 25);

    fmpz_mat_randbits(A, state, 40);
    fmpz_mat_randbits(B, state, 70);

    fmpz_mat_mul(C, A, B);

    fmpz_mat_mul_classical(D, A, B);
    flint_printf("classical:  %d\n", fmpz_mat_equal(C, D));
    fmpz_mat_mul_multi_mod(D, A, B);
    flint_printf("multi_mod:  %d\n", fmpz_mat_equal(C, D));
    fmpz_mat_mul_strassen(D, A, B);
    flint_printf("strassen:   %d\n", fmpz_mat_equal(C, D));
    fmpz_mat_mul_waksman(D, A, B);
    flint_printf("waksman:    %d\n", fmpz_mat_equal(C, D));

    /* these may decline; then D is left for fmpz_mat_mul */
    ok = fmpz_mat_mul_blas(D, A, B);
    flint_printf("blas:       %d\n", !ok || fmpz_mat_equal(C, D));
    ok = fmpz_mat_mul_fft_small(D, A, B);
    flint_printf("fft_small:  %d\n", !ok || fmpz_mat_equal(C, D));

    flint_printf("max bits of C: %wd\n", fmpz_mat_max_bits(C));

    /* square a 2 x 2 matrix in place: Fibonacci numbers */
    fmpz_mat_clear(A);
    fmpz_mat_init(A, 2, 2);
    fmpz_one(fmpz_mat_entry(A, 0, 0));
    fmpz_one(fmpz_mat_entry(A, 0, 1));
    fmpz_one(fmpz_mat_entry(A, 1, 0));
    fmpz_mat_pow(A, A, 100);
    flint_printf("F(100) = %{fmpz}\n", fmpz_mat_entry(A, 0, 1));

    fmpz_mat_clear(A);
    fmpz_mat_clear(B);
    fmpz_mat_clear(C);
    fmpz_mat_clear(D);
    flint_rand_clear(state);
    return 0;
}
```
```output
classical:  1
multi_mod:  1
strassen:   1
waksman:    1
blas:       1
fft_small:  1
max bits of C: -114
F(100) = 354224848179261915075
```
