---
title: Multiplication modulo 2^N + 1
functions: fft_mulmod_2expp1 fft_adjust_limbs fft_adjust_limbs_sd_fft FLINT_FFT_SD_FFT_MIN_RING _fft_mulmod_2expp1 fft_naive_convolution_1
---

The pointwise products of a Schönhage–Strassen convolution are products
modulo $2^N + 1$. `fft_mulmod_2expp1` computes them, by a full product
followed by a reduction for $N$ up to 8192 bits and by a negacyclic
convolution, the original idea of Schönhage and Strassen [@SS1971], for
larger $N$. Since the negacyclic convolution needs $N$ to have a special
form, the outer convolution must choose its coefficient ring with
`fft_adjust_limbs` (or `fft_adjust_limbs_sd_fft`, which also considers
the faster pointwise engine of `fft_small`).

## Usage {#usage}

The operands are reduced residues of $\mathit{limbs} + 1$ limbs, where
$\mathit{limbs} = nw/\bits$: only the product $nw$ matters, and
`fft_convolution` passes the transform parameters $n$ and $w$ of the
outer transform. The function accepts any $N$ up to
`FFT_MULMOD_2EXPP1_CUTOFF` limbs, which is 128 on all targets. Beyond
that, $N$ must be a value produced by `fft_adjust_limbs`; for other
values the result is silently wrong. For example, `fft_adjust_limbs`
maps 129 to 144, 200 to 208, 1000 to 1024 and 12345 to 12544, and
products modulo $2^{64 \cdot 129} + 1$ or $2^{64 \cdot 1000} + 1$ come
out wrong. Every power of two and every multiple of three by a power of
two above 128 is returned unchanged.

`tt` must provide $2 \cdot \mathit{limbs}$ limbs of scratch space for the
full product of the basecase; FLINT's callers pass
$2(\mathit{limbs} + 1)$ limbs. The output may alias either input.

## Algorithm {#algorithm}

If one of the inputs is $2^N \equiv -1$, the result is the negation of
the other. Up to the cutoff, `flint_mpn_mulmod_2expp1_basecase` computes
the full $2N$-bit product with `flint_mpn_mul_n` (or `flint_mpn_sqr`) and
subtracts the high half from the low half. Beyond the cutoff,
`_fft_mulmod_2expp1` cuts both operands into $2n$ pieces of
$b = N/(2n)$ bits, with $n = 2^{\mathit{depth}}$ and
$\mathit{depth} = \lfloor \lceil \log_2 N \rceil / 2 \rfloor - k$, where
$k$ comes from the tuning table `MULMOD_TAB` (between 1 and 4). Since
$2^{2nb} = 2^N \equiv -1$, the product modulo $2^N + 1$ is the
negacyclic convolution of the pieces, evaluated at $2^b$. The
convolution is computed with `fft_negacyclic` and `ifft_negacyclic` in
the ring $\Z/(2^{2b} + 1)\Z$, with `flint_mpn_mulmod_2expp1_basecase` for
its pointwise products, and, to recover coefficients that exceed $2^{2b}$,
combined by the Chinese remainder theorem with the same convolution
modulo $B$, computed by `fft_naive_convolution_1` from the lowest limb of
each piece. The proof that this gives the exact coefficients is in
[Multiplication modulo $2^N + 1$](algorithms/fft/schonhage-strassen#mulmod).
Some values:

| $N$ (bits) | pieces $2n$ | piece bits $b$ | inner ring bits $2b$ |
|---:|---:|---:|---:|
| 16384 | 16 | 1024 | 2048 |
| 65536 | 32 | 2048 | 4096 |
| 262144 | 128 | 2048 | 4096 |
| 1048576 | 256 | 4096 | 8192 |
| 4194304 | 512 | 8192 | 16384 |

(generic tuning; the inner products are handled by the basecase, so
there is a single level of recursion).

`fft_adjust_limbs` rounds `limbs` up so that the depth $d_1$ chosen for
`limbs` and for the next power of two divides it correctly: the limb
count becomes a multiple of $2^{d_1 + 1}$ and the bit count a multiple of
$2^{2 d_1}$, where $d_1$ is the larger of the two depths.

## The negacyclic engine of fft_small {#sd-fft}

In builds with `fft_small`, `sd_fft_mpn_mulmod_2expp1` computes the same
products with a negacyclic convolution modulo three to eight word-size
primes in floating point, for rings $N' = mb$ with $m \ge 16$ a power of
two and $b$ one of 18 instantiated digit sizes between 64 and 192 bits.
`fft_adjust_limbs_sd_fft(tight, n_outer, m_out)` looks for such a ring
for a convolution with transform parameter $n$ = `n_outer` and
coefficients of at least `tight` limbs:

- $m$ is the least power of two $\ge 16$ with $mb \ge \bits \cdot \mathit{tight}$;
- $mb$ must be a multiple of $\bits$, at least
  `FLINT_FFT_SD_FFT_MIN_RING` $= 16384$ bits, and a multiple of
  $2n$, so that the outer transform has an even number $w$ of bits
  per point;
- the number of primes needed for $2b + \log_2 m + 1$ bits must be
  available, and the ring may exceed `fft_adjust_limbs(tight)` by at most
  35, 25, 15, 12, 5 or 2 percent for 3, 4, …, 8 primes;
- among the admissible rings, the one with the least estimated cost
  $\mathit{np} \cdot m \cdot (\log_2 m + 4 + b/16)$ is chosen.

It returns the limb count $mb/\bits$ and stores $m$ in `*m_out`, or
returns `fft_adjust_limbs(tight)` with `*m_out` $= 0$ if no ring
qualifies or `fft_small` is not available. The caps and weights are
tuning parameters. `fmpz_poly` and `fmpzi` select their rings this way
and create a context with `sd_fft_mpn_mulmod_2expp1_ctx_init` when
$m \ne 0$; see the second example in [](fft/convolution#examples).

## Example {#example}

A product modulo $2^N + 1$ with $N$ near 64000 bits, checked against
GMP.

```c example name=fft_mulmod_2expp1
#include <gmp.h>
#include <flint/flint.h>
#include <flint/fft.h>
#include <flint/mpn_extras.h>

int main(void)
{
    /* multiply modulo 2^N + 1 with N a multiple of 64 near 64000 bits */
    slong limbs = fft_adjust_limbs(1000);
    mp_limb_t * a, * b, * r, * tt;
    mpz_t A, B, R, P;

    flint_printf("limbs: 1000 -> %wd\n", limbs);

    a = flint_calloc(limbs + 1, sizeof(mp_limb_t));
    b = flint_calloc(limbs + 1, sizeof(mp_limb_t));
    r = flint_calloc(limbs + 1, sizeof(mp_limb_t));
    tt = flint_malloc(2 * (limbs + 1) * sizeof(mp_limb_t));

    mpz_init(A); mpz_init(B); mpz_init(R); mpz_init(P);
    mpz_ui_pow_ui(P, 2, limbs * FLINT_BITS);
    mpz_add_ui(P, P, 1);
    mpz_ui_pow_ui(A, 3, 40000);
    mpz_ui_pow_ui(B, 5, 27000);
    mpz_mod(A, A, P);
    mpz_mod(B, B, P);

    /* reduced operands: limbs + 1 limbs, top limb 0 (or 1 for 2^N) */
    mpz_export(a, NULL, -1, sizeof(mp_limb_t), 0, 0, A);
    mpz_export(b, NULL, -1, sizeof(mp_limb_t), 0, 0, B);

    /* n * w = N; here n = limbs and w = FLINT_BITS */
    fft_mulmod_2expp1(r, a, b, limbs, FLINT_BITS, tt);

    fermat_to_mpz(R, r, limbs);
    mpz_mul(A, A, B);
    mpz_mod(A, A, P);
    flint_printf("equal to mpz: %d\n", mpz_cmp(A, R) == 0);

    flint_free(a); flint_free(b); flint_free(r); flint_free(tt);
    mpz_clear(A); mpz_clear(B); mpz_clear(R); mpz_clear(P);
    return 0;
}
```
```output
limbs: 1000 -> 1024
equal to mpz: 1
```
