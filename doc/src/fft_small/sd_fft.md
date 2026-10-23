---
title: Transforms modulo one prime
functions: sd_fft_trunc sd_ifft_trunc sd_fft_trunc_tr sd_fft_ctx_point_mul
---

`sd_fft_trunc` and `sd_ifft_trunc` are the number-theoretic transforms
that all of `fft_small` is built on. They work in place on an array of
doubles holding integers, modulo a prime $p < 2^{50}$ with $2^L \mid
p - 1$, and are truncated: only part of the input is nonzero and only
part of the output is needed. Using them directly makes sense for
convolutions modulo one of the context primes, or modulo another suitable
prime with a context of its own.

## What they compute {#semantics}

Let $w_0, w_1, w_2, \ldots$ be the roots of unity of `sd_fft_ctx_w`: the
first $2^L$ of them are the $2^L$-th roots of unity in bit-reversed order,
$w_i = \omega^{\operatorname{rev}_L(i)}$ for a primitive $2^L$-th root
$\omega$. For the polynomial $f = \sum_{j < \mathtt{itrunc}} d_j x^j$
formed by the first `itrunc` entries, `sd_fft_trunc` stores a value
congruent to $f(w_i)$ at position `sd_fft_ctx_trunc_index(L, i)`, for
$i < \mathtt{otrunc}$. That position is $i$ for $L \le 3$; for $L \ge 4$
each group of 16 values is stored as a transposed $4 \times 4$ matrix,
because the innermost kernel saves its final transposition. Pointwise
operations do not care about the order, and the inverse transform expects
the same one.

`sd_ifft_trunc` inverts the transform up to the factor $2^L$: if the first
`trunc` positions hold the values of the forward transform of a
polynomial $g$ of length at most `trunc`, it replaces them by values
congruent to $2^L g_0, \ldots, 2^L g_{\mathtt{trunc}-1}$. The truncated
inverse needs no other values, because the missing evaluations are
determined by the known zero coefficients [@vdH2004].

`sd_fft_ctx_point_mul` multiplies the transforms pointwise and by a
constant $m$ in one pass. Taking $m = 2^{-L} \bmod p$ makes the inverse
transform return the product itself.

`sd_fft_trunc_tr` and `sd_ifft_trunc_tr` are the transposes of the two
transforms, with the same parameters: `sd_fft_trunc_tr` reads `otrunc`
values on the evaluation side, in the order above, and writes `itrunc`
coefficients $\sum_{i < \mathtt{otrunc}} v_i w_i^t$, $t < \mathtt{itrunc}$,
where $v_i$ is the entry at position `sd_fft_ctx_trunc_index(L, i)`. See
[](#transposed).

## Requirements {#requirements}

- `d` must be aligned to 32 bytes and hold $2^L$ doubles.
- Entries are integers. The forward transform accepts and returns entries
  in $(-3p, 3p)$. The inverse transform accepts entries in $(-2p, 2p)$,
  which `sd_fft_ctx_point_mul` produces from forward outputs, and returns
  entries in $(-3p, 3p)$. Reduce with `nmod_set_si` or the vector
  reductions of `machine_vectors.h` to obtain residues.
- Truncations are rounded: for $L \le 8$ the forward transform zeroes the
  entries from `itrunc` on and computes all $2^L$ values, and the inverse
  transform ignores `trunc` and needs all $2^L$ values. For $L > 8$ both
  work on blocks of 256 entries: `itrunc`, `otrunc` and `trunc` are
  rounded up to multiples of 256, the forward transform zeroes the
  entries from `itrunc` to the next multiple of 256, and entries beyond
  the rounded output truncation are neither read nor written.
- `sd_fft_ctx_point_mul` touches all $2^L$ entries of both arrays, so
  entries beyond the truncation must hold some finite values; it requires
  $L \ge 4$ and $0 \le m \le p$.
- The transposed transforms have the rounding of the transforms they
  transpose. `sd_ifft_trunc_tr` accepts and returns entries in
  $(-3p, 3p)$, like `sd_fft_trunc`; `sd_fft_trunc_tr` accepts entries in
  $(-2p, 2p)$, like `sd_ifft_trunc`, and returns entries in $(-2p, 2p)$.
- The context must cover depth $L$, which all four transforms ensure by
  calling `sd_fft_ctx_fit_depth`, and $2^L$ must divide $p - 1$.

## Algorithm {#algorithm}

The transform is a fast Fourier transform [@CT1965] over $\Z/p\Z$
[@Pol1971], computed as a tree of reductions of $f$ modulo
$x^{2^k} - w_j$ [@Ber2008]: the butterfly $(a, b) \mapsto (a + sb,
a - sb)$ splits a residue modulo $x^{2h} - s^2$ into its residues modulo
$x^h - s$ and $x^h + s$. Two levels are merged into radix-4 butterflies. Above 256
points it is organised as Bailey's matrix algorithm [@Bai1990]: the
length $2^k$ blocks are viewed as a matrix of $2^{k_1}$ rows of
$2^{k_2}$ blocks, $k_1 = \lfloor k/2 \rfloor$, and column transforms on
strided blocks are followed by row transforms with twisted roots, each
recursively. Within a block, fixed-size kernels keep the data in
registers. The inverse transform runs the same structure backwards with
Gentleman–Sande butterflies [@GS1966]. Each butterfly uses the
floating-point modular arithmetic of [](algorithms/simd/fp-modular-arithmetic),
with lazy reductions whose ranges are proved on the algorithm page.

::: summary algorithms/fft/small-prime-fft :::

## Transposed transforms {#transposed}

A transform computed by a linear program of butterflies has a transpose
computed by the same butterflies in reverse order, each replaced by its
transpose, at exactly the same cost (the transposition principle,
[](algorithms/polynomials/middle-product#tellegen)). The transpose of the
butterfly $(a, b) \mapsto (a + wb, a - wb)$ of the forward transform is
$(a, b) \mapsto (a + b, w(a - b))$, the shape of the inverse butterflies
with the forward twiddles, and conversely; truncation transposes into
truncation, with the roles of `itrunc` and `otrunc` exchanged.

Their use is the transposed convolution. For a fixed $a$ of length $a_n$,
multiplication by $a$ of polynomials of length $z_n$ is, for a transform
truncation $T \ge a_n + z_n - 1$, the composition
$b \mapsto I(m\,\hat a \cdot F(b))$ of the forward transform $F$ of $b$
(`itrunc` $= z_n$, `otrunc` $= T$), the pointwise product with
$m\,\hat a$, and the inverse transform $I$ with truncation $T$, followed by
taking the first $a_n + z_n - 1$ coefficients. Its transpose is
$v \mapsto F^T(m\,\hat a \cdot I^T(v))$ for $v$ of length
$a_n + z_n - 1$ padded with zeros to $T$, and by
[](algorithms/polynomials/middle-product#lem-transpose) it computes
$z_j = \sum_i a_i v_{i+j}$ for $j < z_n$: a middle product, at the cost of
a product of lengths $a_n$ and $z_n$.
`_nmod_poly_mul_tr_mpn_ctx` is this computation modulo several primes;
see [](algorithms/fft/small-prime-fft#transposed).

## Examples {#examples}

The positions of the first values of a transform of length 16, and the
inverse transform:

```c example name=fft_small_sd_fft_order requires=fft_small
#include <flint/flint.h>
#include <flint/nmod.h>
#include <flint/fft_small.h>

int main(void)
{
    mpn_ctx_struct * R = get_default_mpn_ctx();
    sd_fft_ctx_struct * Q = R->ffts + 0;      /* p = 63*2^44 + 1 */
    nmod_t mod = Q->mod;
    ulong L = 4, N = 16, i, j, k, fx, x, inv16;
    ulong f[3] = { 1, 2, 3 };                  /* f = 1 + 2x + 3x^2 */
    double * d = flint_aligned_alloc(32, N * sizeof(double));

    for (i = 0; i < 3; i++)
        d[i] = f[i];

    /* forward transform: f(w_i) lands at position trunc_index(L, i) */
    sd_fft_trunc(Q, d, L, 3, N);

    for (i = 0; i < N; i++)
    {
        x = nmod_set_si((slong) sd_fft_ctx_w(Q, i), mod);
        fx = 0;
        for (j = 3; j-- > 0; )
            fx = nmod_add(nmod_mul(fx, x, mod), f[j], mod);
        k = sd_fft_ctx_trunc_index(L, i);
        if (nmod_set_si((slong) d[k], mod) != fx)
            flint_printf("mismatch at i = %wu\n", i);
        if (i < 6)
            flint_printf("w_%wu: position %2wu\n", i, k);
    }

    /* inverse transform: returns 16 f, to be divided by 2^L */
    sd_ifft_trunc(Q, d, L, N);
    inv16 = nmod_inv(16, mod);
    for (i = 0; i < 4; i++)
        flint_printf("%wu ", nmod_mul(nmod_set_si((slong) d[i], mod), inv16, mod));
    flint_printf("\n");

    flint_aligned_free(d);
    flint_cleanup_master();
    return 0;
}
```
```output
w_0: position  0
w_1: position  4
w_2: position  8
w_3: position 12
w_4: position  1
w_5: position  5
1 2 3 0
```

A product of polynomials of lengths 300 and 200 modulo the prime with
truncated transforms of length $2^{10}$, of which only the first 512
values are computed:

```c example name=fft_small_sd_fft_mul requires=fft_small
#include <flint/flint.h>
#include <flint/nmod.h>
#include <flint/nmod_vec.h>
#include <flint/nmod_poly.h>
#include <flint/fft_small.h>

int main(void)
{
    sd_fft_ctx_struct * Q = get_default_mpn_ctx()->ffts + 0;
    nmod_t mod = Q->mod;
    slong an = 300, bn = 200, zn = an + bn - 1, i;
    ulong L = 10, m;
    nn_ptr a = _nmod_vec_init(an), b = _nmod_vec_init(bn);
    nn_ptr c = _nmod_vec_init(zn), z = _nmod_vec_init(zn);
    double * x = flint_aligned_alloc(32, 1024 * sizeof(double));
    double * y = flint_aligned_alloc(32, 1024 * sizeof(double));

    for (i = 0; i < an; i++)
        a[i] = nmod_pow_ui(i + 1, 7, mod);
    for (i = 0; i < bn; i++)
        b[i] = nmod_pow_ui(i + 2, 9, mod);

    /* residues in [0, p) are valid inputs; the unused tail is zeroed
       because sd_fft_ctx_point_mul reads all 2^L entries */
    for (i = 0; i < 1024; i++)
        x[i] = y[i] = 0.0;
    for (i = 0; i < an; i++)
        x[i] = (double) a[i];
    for (i = 0; i < bn; i++)
        y[i] = (double) b[i];

    /* only 512 of the 1024 values are needed, since zn <= 512 */
    sd_fft_trunc(Q, x, L, an, 512);
    sd_fft_trunc(Q, y, L, bn, 512);

    /* x = x * y / 2^L pointwise, then the inverse transform */
    m = nmod_inv(nmod_pow_ui(2, L, mod), mod);
    sd_fft_ctx_point_mul(Q, x, y, m, L);
    sd_ifft_trunc(Q, x, L, 512);

    for (i = 0; i < zn; i++)
        z[i] = nmod_set_si((slong) x[i], mod);

    _nmod_poly_mul(c, a, an, b, bn, mod);
    flint_printf("product modulo p: %s\n",
                 _nmod_vec_equal(z, c, zn) ? "agrees" : "WRONG");

    _nmod_vec_clear(a);
    _nmod_vec_clear(b);
    _nmod_vec_clear(c);
    _nmod_vec_clear(z);
    flint_aligned_free(x);
    flint_aligned_free(y);
    flint_cleanup_master();
    return 0;
}
```
```output
product modulo p: agrees
```

The transposed product $z_j = \sum_i a_i v_{i+j}$, $j < 400$, of $a$ of
length 300 and $v$ of length 699 modulo the prime, with transforms of
length $2^{10}$ truncated to 768 points:

```c example name=fft_small_sd_fft_transposed requires=fft_small
#include <flint/flint.h>
#include <flint/nmod.h>
#include <flint/nmod_vec.h>
#include <flint/fft_small.h>

int main(void)
{
    sd_fft_ctx_struct * Q = get_default_mpn_ctx()->ffts + 0;
    nmod_t mod = Q->mod;
    slong an = 300, zn = 400, vn = an + zn - 1, i, j;
    ulong L = 10, T = 768, m, s;
    nn_ptr a = _nmod_vec_init(an), v = _nmod_vec_init(vn);
    nn_ptr z = _nmod_vec_init(zn);
    double * x = flint_aligned_alloc(32, 1024 * sizeof(double));
    double * y = flint_aligned_alloc(32, 1024 * sizeof(double));
    int ok = 1;

    for (i = 0; i < an; i++)
        a[i] = nmod_pow_ui(i + 1, 7, mod);
    for (i = 0; i < vn; i++)
        v[i] = nmod_pow_ui(i + 2, 9, mod);

    for (i = 0; i < 1024; i++)
        x[i] = y[i] = 0.0;
    for (i = 0; i < an; i++)
        x[i] = (double) a[i];
    for (i = 0; i < vn; i++)
        y[i] = (double) v[i];

    /* a is transformed as for a product of length T */
    sd_fft_trunc(Q, x, L, an, T);

    /* the transposed pipeline: ifft_tr, pointwise product, fft_tr */
    m = nmod_inv(nmod_pow_ui(2, L, mod), mod);
    sd_ifft_trunc_tr(Q, y, L, T);
    sd_fft_ctx_point_mul(Q, y, x, m, L);
    sd_fft_trunc_tr(Q, y, L, zn, T);

    for (j = 0; j < zn; j++)
    {
        s = 0;
        for (i = 0; i < an; i++)
            s = nmod_addmul(s, a[i], v[i + j], mod);
        z[j] = nmod_set_si((slong) y[j], mod);
        ok = ok && (z[j] == s);
    }
    flint_printf("transposed product: %s\n", ok ? "agrees" : "WRONG");

    _nmod_vec_clear(a);
    _nmod_vec_clear(v);
    _nmod_vec_clear(z);
    flint_aligned_free(x);
    flint_aligned_free(y);
    flint_cleanup_master();
    return 0;
}
```
```output
transposed product: agrees
```
