---
title: Convolutions of large integers
functions: fft_convolution fft_convolution_sd_fft fft_convolution_choose_mid
---

`fft_convolution` multiplies two polynomials whose coefficients are
residues modulo $p = 2^N + 1$, $N = \mathit{limbs} \cdot \bits$. When
the coefficients of the true integer product are known to be smaller
than $p$ (or than $p/2$ in absolute value), it computes the product of
two polynomials with integer coefficients exactly. This is how `fmpz_poly`
implements Schönhage–Strassen multiplication, and with integer
coefficients obtained by cutting integers into chunks, it is also the
core of `mul_mfa_truncate_sqrt2`. The guide [](fft/using) explains the
whole procedure.

## Parameters and layout {#layout}

The transform has length $4n$, $n = 2^{\mathit{depth}}$, and only
`trunc` of its points are computed. To multiply polynomials of lengths
$\ell_1$ and $\ell_2$, take $\mathit{trunc} = \ell_1 + \ell_2 - 1$ and
$\mathit{depth} = \lceil \log_2 \mathit{trunc} \rceil - 2$, so that
$2n < \mathit{trunc} \le 4n$ as required (for $\mathit{trunc} \ge 3$).
The ring size $N$ must be a multiple of $n$, so that the transform has
$w = N/n \ge 1$ bits per point, and if $\mathit{limbs} > 128$ it must
be a value returned by `fft_adjust_limbs` (or `fft_adjust_limbs_sd_fft`)
because the pointwise products use `fft_mulmod_2expp1`.

- `ii` and `jj` are arrays of $4n$ pointers to blocks of
  $\mathit{limbs} + 1$ limbs. The first $\ell_1$, respectively $\ell_2$,
  entries hold the input coefficients as $(\mathit{limbs} + 1)$-limb
  two's complement integers (negative values are allowed); all other
  entries must be zero. Only the entries below `trunc`, rounded up to an
  even number or, for $\mathit{depth} > 6$, to a multiple of
  $2^{\lfloor \mathit{depth}/2 \rfloor + 1}$, are read, but zeroing all
  of them is simplest.
- On return, `ii[k]` holds the coefficient $c_k$ of the product, reduced
  to $[0, p)$, for $k < \mathit{trunc}$. If the true coefficients satisfy
  $|c_k| < p/2$, then $c_k$ is the residue if it is less than $p/2$ and
  the residue minus $p$ otherwise.
- `jj` is overwritten by its transform. Passing the same array twice
  (`ii == jj`) computes a square with one forward transform less.
- `t1`, `t2` and `s1` are arrays of at least `flint_get_num_threads()`
  pointers to $\mathit{limbs} + 1$ limbs each, and `tt` an array of as
  many pointers to $2(\mathit{limbs} + 1)$ limbs each. For
  $\mathit{depth} \le 6$ only the first entries are used.

The bound on $N$ follows from the sizes of the inputs: if
$|a_i| < 2^{b_1}$ and $|b_j| < 2^{b_2}$, then
$|c_k| < \min(\ell_1, \ell_2) \, 2^{b_1 + b_2}$, so

$$
N \ge b_1 + b_2 + \lceil \log_2 \min(\ell_1, \ell_2) \rceil + s
$$

suffices, with $s = 1$ if any coefficient is negative and $s = 0$
otherwise. `_fmpz_poly_mulmid_SS` uses exactly this bound, rounds it up
to a multiple of $n$ and of $\bits$, and adjusts the limb count with
`fft_adjust_limbs_sd_fft`.

If the product has more than `trunc` coefficients, the result is
meaningless, except for $\mathit{trunc} = 4n$, where `fft_convolution`
computes the cyclic product modulo $x^{4n} - 1$ (each coefficient $c_k$
with $k \ge 4n$ is added to $c_{k - 4n}$). `fft_convolution_choose_mid`
exploits this.

## Algorithm {#algorithm}

For $\mathit{depth} \le 6$, `fft_convolution` uses `fft_truncate_sqrt2`
on both operands, multiplies the transformed coefficients with
`fft_mulmod_2expp1`, applies `ifft_truncate_sqrt2` and divides by $4n$.
For larger depths it uses the [matrix Fourier algorithm](algorithms/fft/matrix-fourier)
with $n_1 = 2^{\lfloor \mathit{depth}/2 \rfloor}$ columns, in three
threaded stages: `fft_mfa_truncate_sqrt2_outer` (sqrt2 layer and column
transforms) on each operand, then `fft_mfa_truncate_sqrt2_inner`, which
handles one row at a time (row transforms of both operands, pointwise
products, inverse row transform) while the row is in cache, and finally
`ifft_mfa_truncate_sqrt2_outer`, which also divides by $4n$. The total
cost is that of three [truncated transforms](algorithms/fft/truncated-fft)
of `trunc` points plus `trunc` multiplications modulo $2^N + 1$, each a
[Schönhage–Strassen product](algorithms/fft/schonhage-strassen#mulmod)
of its own for large $N$.

`fft_convolution_sd_fft` computes the pointwise products with
`sd_fft_mpn_mulmod_2expp1` from `fft_small` instead, a negacyclic
convolution modulo several word-size primes
([](fft_small/mulmod_2expp1)), which is about twice as
fast on a recent x86-64 processor for rings of $16384$ bits and more. Its context must be initialised with
`sd_fft_mpn_mulmod_2expp1_ctx_init` for $N = \mathit{limbs} \cdot \bits$
and the digit count $m$ that `fft_adjust_limbs_sd_fft` returned together
with `limbs`. For large depths it keeps the fused inner stage only when
there are at least eight work units (columns plus truncated rows) per
thread; otherwise it performs the transforms in separate single-threaded
passes and only the pointwise products, which dominate, in threads.

## Partial products {#mid}

When only the coefficients with indices in $[\mathit{nlo}, \mathit{nhi})$
of a product of length $\ell = \ell_1 + \ell_2 - 1$ are wanted, a
cyclic transform can be shorter than the truncated one. In the product
modulo $x^M - 1$, coefficient $k$ receives $c_k + c_{k+M}$. If
$M \ge \mathit{nhi}$ and $M \ge \ell - \mathit{nlo}$, every wrapped
term $c_{k+M}$ with $k \ge \mathit{nlo}$ has $k + M \ge \ell$ and is zero,
so the window is unaffected. `fft_convolution_choose_mid` takes the
least power of two $M = 4 \cdot 2^{\mathit{depth}}$ with this property
and returns it as `trunc` if it is smaller than the number of points the
truncated transform of length $\ell$ would process; otherwise it returns
$\mathit{depth} = \lceil \log_2 \ell \rceil - 2$ (at least 0) with
$\mathit{trunc} = \ell$. For the middle product of $2m$ by $m$
coefficients, where the wanted window has length $m$, this gives $M = 2m$
points instead of about $3m$ when $2m$ is a power of two.

## Examples {#examples}

The product of two polynomials of length 300 with signed 31-bit
coefficients. The coefficients of the product are below
$300 \cdot 2^{60} < 2^{69}$ in absolute value, so $N = 128$ would do,
but $N$ must be a multiple of $n = 256$, so the example uses four limbs.

```c example name=fft_convolution_signed
#include <gmp.h>
#include <flint/flint.h>
#include <flint/fft.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/mpn_extras.h>

/* one coefficient c with |c| < 2^63, as a (limbs + 1)-limb two's complement value */
static void set_si(mp_limb_t * x, slong c, slong limbs)
{
    flint_mpn_store(x, limbs + 1, c < 0 ? ~UWORD(0) : 0);
    x[0] = (mp_limb_t) c;
}

/* a reduced residue modulo 2^N + 1, read in the symmetric range */
static void get_fmpz(fmpz_t c, mp_limb_t * x, slong limbs, const mpz_t p)
{
    mpz_t t;
    mpz_init(t);
    fermat_to_mpz(t, x, limbs);
    if (mpz_sizeinbase(t, 2) >= limbs * FLINT_BITS)   /* t >= 2^(N-1) */
        mpz_sub(t, t, p);
    fmpz_set_mpz(c, t);
    mpz_clear(t);
}

int main(void)
{
    slong len = 300, len_out = 2 * len - 1;
    slong depth = FLINT_CLOG2(len_out) - 2;     /* 4n >= len_out > 2n */
    slong n = WORD(1) << depth;                  /* 256 */
    slong limbs = 4, size = limbs + 1;           /* N = 256 >= 70 bits, n | N */
    slong nt = flint_get_num_threads(), i;
    mp_limb_t ** ii, ** jj, ** t1, ** t2, ** s1, ** tt, * mem;
    fmpz_poly_t f, g, h, h2;
    mpz_t p;

    fmpz_poly_init(f); fmpz_poly_init(g); fmpz_poly_init(h); fmpz_poly_init(h2);
    for (i = 0; i < len; i++)
    {
        fmpz_poly_set_coeff_si(f, i, (i * 1103515245 + 12345) % (WORD(1) << 30) - (WORD(1) << 29) * (i % 3));
        fmpz_poly_set_coeff_si(g, i, (i * 69069 + 1) % (WORD(1) << 30) - (WORD(1) << 30) * (i % 2));
    }

    ii = flint_malloc((8 * n + 4 * nt) * sizeof(mp_limb_t *));
    jj = ii + 4 * n;
    t1 = jj + 4 * n; t2 = t1 + nt; s1 = t2 + nt; tt = s1 + nt;
    mem = flint_malloc((8 * n + 5 * nt) * size * sizeof(mp_limb_t));
    for (i = 0; i < 8 * n; i++)
        ii[i] = mem + i * size;
    for (i = 0; i < nt; i++)
    {
        t1[i] = mem + (8 * n + 5 * i) * size;
        t2[i] = t1[i] + size;
        s1[i] = t2[i] + size;
        tt[i] = s1[i] + size;               /* 2 (limbs + 1) limbs */
    }

    for (i = 0; i < 4 * n; i++)
    {
        set_si(ii[i], i < len ? fmpz_get_si(f->coeffs + i) : 0, limbs);
        set_si(jj[i], i < len ? fmpz_get_si(g->coeffs + i) : 0, limbs);
    }

    fft_convolution(ii, jj, depth, limbs, len_out, t1, t2, s1, tt);

    mpz_init(p);
    mpz_ui_pow_ui(p, 2, limbs * FLINT_BITS);
    mpz_add_ui(p, p, 1);
    fmpz_poly_fit_length(h, len_out);
    for (i = 0; i < len_out; i++)
        get_fmpz(h->coeffs + i, ii[i], limbs, p);
    _fmpz_poly_set_length(h, len_out);
    _fmpz_poly_normalise(h);

    fmpz_poly_mul(h2, f, g);
    flint_printf("depth = %wd, trunc = %wd, equal: %d\n", depth, len_out, fmpz_poly_equal(h, h2));
    flint_printf("coefficient 299: ");
    fmpz_print(h->coeffs + 299);
    flint_printf("\n");

    mpz_clear(p);
    flint_free(mem); flint_free(ii);
    fmpz_poly_clear(f); fmpz_poly_clear(g); fmpz_poly_clear(h); fmpz_poly_clear(h2);
    return 0;
}
```
```output
depth = 8, trunc = 599, equal: 1
coefficient 299: -889640045340373524
```

A product with 8000-bit coefficients, where the ring has 16384 bits and
the pointwise products can use the negacyclic engine of `fft_small`. Here
`fft_adjust_limbs_sd_fft` rounds 251 limbs up to 256 and returns
$m = 256$ digits of 64 bits. In builds without `fft_small` it returns
$m = 0$ and the program uses `fft_convolution`.

```c example name=fft_convolution_sd_fft
#include <gmp.h>
#include <flint/flint.h>
#include <flint/fft.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/mpn_extras.h>
#if FLINT_HAVE_FFT_SMALL
#include <flint/fft_small.h>
#endif

int main(void)
{
    slong len = 40, len_out = 2 * len - 1, i, m;
    slong depth = FLINT_CLOG2(len_out) - 2, n = WORD(1) << depth;   /* n = 32 */
    slong bits = 8000, res_bits, limbs, size;
    slong nt = flint_get_num_threads();
    mp_limb_t ** ii, ** jj, ** t1, ** t2, ** s1, ** tt, * mem;
    fmpz_poly_t f, g, h, h2;
    flint_rand_t state;

    flint_rand_init(state);
    fmpz_poly_init(f); fmpz_poly_init(g); fmpz_poly_init(h); fmpz_poly_init(h2);
    fmpz_poly_randtest_unsigned(f, state, len, bits);
    fmpz_poly_randtest_unsigned(g, state, len, bits);
    fmpz_poly_set_coeff_ui(f, len - 1, 1);
    fmpz_poly_set_coeff_ui(g, len - 1, 1);

    /* product coefficients < len 2^(2 bits): round up to a multiple of n */
    res_bits = 2 * bits + FLINT_CLOG2(len);
    res_bits = ((res_bits + n - 1) / n) * n;
    limbs = fft_adjust_limbs_sd_fft((res_bits - 1) / FLINT_BITS + 1, n, &m);
    size = limbs + 1;

    ii = flint_malloc((8 * n + 4 * nt) * sizeof(mp_limb_t *));
    jj = ii + 4 * n;
    t1 = jj + 4 * n; t2 = t1 + nt; s1 = t2 + nt; tt = s1 + nt;
    mem = flint_calloc((8 * n + 5 * nt) * size, sizeof(mp_limb_t));
    for (i = 0; i < 8 * n; i++)
        ii[i] = mem + i * size;
    for (i = 0; i < nt; i++)
    {
        t1[i] = mem + (8 * n + 5 * i) * size;
        t2[i] = t1[i] + size; s1[i] = t2[i] + size; tt[i] = s1[i] + size;
    }

    for (i = 0; i < len; i++)
    {
        fmpz_get_ui_array(ii[i], limbs, f->coeffs + i);
        fmpz_get_ui_array(jj[i], limbs, g->coeffs + i);
    }

#if FLINT_HAVE_FFT_SMALL
    if (m != 0)
    {
        sd_fft_mpn_mulmod_2expp1_ctx_t ctx;
        sd_fft_mpn_mulmod_2expp1_ctx_init(ctx, get_default_mpn_ctx(),
                                          limbs * FLINT_BITS, m);
        fft_convolution_sd_fft(ctx, ii, jj, depth, limbs, len_out, t1, t2, s1, tt);
        sd_fft_mpn_mulmod_2expp1_ctx_clear(ctx);
    }
    else
#endif
        fft_convolution(ii, jj, depth, limbs, len_out, t1, t2, s1, tt);

    fmpz_poly_fit_length(h, len_out);
    for (i = 0; i < len_out; i++)
        fmpz_set_ui_array(h->coeffs + i, ii[i], limbs);
    _fmpz_poly_set_length(h, len_out);
    fmpz_poly_mul(h2, f, g);
    flint_printf("equal: %d\n", fmpz_poly_equal(h, h2));

    flint_free(mem); flint_free(ii);
    fmpz_poly_clear(f); fmpz_poly_clear(g); fmpz_poly_clear(h); fmpz_poly_clear(h2);
    flint_rand_clear(state);
    return 0;
}
```
```output
equal: 1
```
