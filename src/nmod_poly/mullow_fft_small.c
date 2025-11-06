/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdint.h>
#include <math.h>
#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

#if FLINT_HAVE_FFT_SMALL

#include "fft_small.h"

/* Repacking for tiny moduli: given polynomials

   A = a[0] + a[1] x + a[2] x^2 + ...
   B = b[0] + b[1] x + b[2] x^2 + ...

   modulo n, we can compute A*B via the half-length product

   A2 = (a[0] + a[1] M) + (a[2] + a[3] M) x + ...
   B2 = (b[0] + b[1] M) + (b[2] + b[3] M) x + ...

   modulo some larger FFT_PRIME if M is large enough that the coefficients
   of A*B do not overflow M and small enough that the coefficients of A2*B2
   do not overflow FFT_PRIME.

   It would be slightly more efficient to combine this repacking with
   the conversion from/to doubles in fft_small, but fft_small is complex
   enough as it is, so we do it here as a preprocessing step instead.
*/

/* Important: this should be one of the default primes used
   by the fft_small default mpn context; otherwise, things will be slow. */
#define FFT_PRIME UWORD(1108307720798209)

/* M will be one of these moduli. M does not need to be a power of two; it
   could be something like M' = floor(FFT_PRIME^(1/3)) instead. However,
   it turns out that M17 usually works when M' does, so currently we just try
   one of the following powers of two. */
#define M16 (UWORD(1) << 16)
#define M17 (UWORD(1) << 17)

/* Only use repacking for 1 <= n <= FFT_SMALL_MAX_REPACK. */
#define FFT_SMALL_MAX_REPACK 23

/* Largest length where the packing became impossible in a test run with
   randtest entries; don't jump into the O(n) norm calculations if we have
   something much larger. This check is not needed for correctness;
   it might furthermore want disabling if we care about often
   multiplying structured with much smaller norms than generic ones. */
static const int repack_limit_tab[] = {
    0,
    264497,  /* 1 */
    256720,  /* 2 */
    77049,  /* 3 */
    36925,  /* 4 */
    21385,  /* 5 */
    14099,  /* 6 */
    9872,  /* 7 */
    7337,  /* 8 */
    5565,  /* 9 */
    4481,  /* 10 */
    3609,  /* 11 */
    3084,  /* 12 */
    2562,  /* 13 */
    2129,  /* 14 */
    1859,  /* 15 */
    1687,  /* 16 */
    1446,  /* 17 */
    1264,  /* 18 */
    1160,  /* 19 */
    1035,  /* 20 */
    942,  /* 21 */
    874,  /* 22 */
    725,  /* 23 */
    662,  /* 24 */
    626,  /* 25 */
    587,  /* 26 */
    562,  /* 27 */
    517,  /* 28 */
    490,  /* 29 */
    426,  /* 30 */
    410,  /* 31 */
    363,  /* 32 */
    342,  /* 33 */
    336,  /* 34 */
    321,  /* 35 */
};

/* Alternative table with cutoffs tuned for more sparse polynomials. */
#if 0
    264497,  /* 1 */
    264497,  /* 2 */
    80976,  /* 3 */
    39194,  /* 4 */
    22473,  /* 5 */
    16039,  /* 6 */
    11338,  /* 7 */
    7709,  /* 8 */
    6647,  /* 9 */
    4945,  /* 10 */
    3904,  /* 11 */
    3574,  /* 12 */
    2770,  /* 13 */
    2690,  /* 14 */
    2213,  /* 15 */
    1841,  /* 16 */
    1607,  /* 17 */
    1474,  /* 18 */
    1339,  /* 19 */
    1228,  /* 20 */
    1805,  /* 21 */
    1562,  /* 22 */
    1517,  /* 23 */
    1432,  /* 24 */
    1252,  /* 25 */
    1138,  /* 26 */
    1085,  /* 27 */
    978,  /* 28 */
    933,  /* 29 */
    858,  /* 30 */
    826,  /* 31 */
    711,  /* 32 */
    698,  /* 33 */
    680,  /* 34 */
    620,  /* 35 */
    587,  /* 36 */
    552,  /* 37 */
    494,  /* 38 */
    474,  /* 39 */
    478,  /* 40 */
    406,  /* 41 */
    414,  /* 42 */
    406,  /* 43 */
    396,  /* 44 */
    360,  /* 45 */
#endif

/* Lemire, Kaser & Kurz method to compute x mod d for 32-bit d using
   a single multiplication.

   The #if-ed version is specialized for half-length (~16-bit moduli).
   This is priori correct for all 16-bit x and d, but we actually need x all
   the way up to 17 bits. Fortunately, it can be checked by brute force to be
   correct for the range of d we need (e.g. 17-bit x and 15-bit d are fine).
   It is easy to find counterexamples by choosing larger x or d,
   e.g. x = 130119, d = 33096 fails.

   We currently disable this version because it is actually benchmarks
   worse on Zen 3 for unknown reasons; perhaps manual SIMD code is needed
   for best performance.
*/

#if 1

FLINT_FORCE_INLINE
uint32_t u32_mod_preinv(uint32_t x, uint32_t d, uint64_t dinv)
{
    uint64_t l = dinv * x;

#if !defined(__GNUC__)
    uint64_t hi, lo;
    umul_ppmm(hi, lo, l, d);
    return hi;
#else
    return ((__uint128_t) l * d ) >> 64;
#endif
}

FLINT_FORCE_INLINE
uint64_t u32_preinv(uint32_t d)
{
    return UINT64_C(0xFFFFFFFFFFFFFFFF) / d + 1;
}

#else

FLINT_FORCE_INLINE
uint32_t u32_mod_preinv(uint32_t x, uint32_t d, uint32_t dinv)
{
    FLINT_ASSERT(x <= (1 << 17) && d <= (1 << 15));

    uint32_t l = dinv * x;
    return ((uint64_t) l * d ) >> 32;
}

FLINT_FORCE_INLINE
uint32_t u32_preinv(uint32_t d)
{
    return UINT32_C(0xFFFFFFFF) / d + 1;
}

#endif

#define UNPACK(m) \
    cy = 0; \
    uint32_t dd = mod.n; \
    uint64_t ddinv = u32_preinv(dd); \
    for (i = 0; 2 * i + 1 < zn; i++) \
    { \
        c = z2[i]; \
        c0 = c % m; \
        c1 = (c / m) % m;  \
        c2 = c / (m * m); \
        d = u32_mod_preinv(c0 + cy, dd, ddinv); \
        z[2 * i] = d; \
        d = u32_mod_preinv(c1, dd, ddinv); \
        z[2 * i + 1] = d; \
        cy = c2; \
    } \
    if (zn % 2) \
    { \
        if (zn > 2 * zn2) \
            c0 = 0; \
        else \
            c0 = z2[zn2 - 1] % M; \
        d = u32_mod_preinv(c0 + cy, dd, ddinv); \
        z[zn - 1] = d; \
    }


static int
_nmod_poly_mullow_fft_small_repack_m(nn_ptr z, nn_srcptr a, slong an, nn_srcptr b, slong bn, slong zn, ulong M, nmod_t mod)
{
    TMP_INIT;
    ulong *a2, *b2, *z2;
    slong an2, bn2, zn2;
    ulong d, cy, c, c0, c1, c2;
    nmod_t mod2;
    slong i;
    int squaring;

    FLINT_ASSERT(M == M16 || M == M17);
    FLINT_ASSERT(mod.n <= FFT_SMALL_MAX_REPACK);

    squaring = (a == b) && (an == bn);

    TMP_START;

    an2 = (an + 1) / 2;
    bn2 = (bn + 1) / 2;
    zn2 = (zn + 1) / 2;
    zn2 = FLINT_MIN(zn2, an2 + bn2 - 1);

    if (squaring)
    {
        a2 = TMP_ALLOC((an2 + zn2) * sizeof(ulong));
        b2 = a2;
        z2 = a2 + an2;
    }
    else
    {
        a2 = TMP_ALLOC((an2 + bn2 + zn2) * sizeof(ulong));
        b2 = a2 + an2;
        z2 = b2 + bn2;
    }

    if (M == M16)
    {
        for (i = 0; i + 1 < an; i += 2)
            a2[i / 2] = a[i] | (a[i + 1] << 16);

        if (!squaring)
            for (i = 0; i + 1 < bn; i += 2)
                b2[i / 2] = b[i] | (b[i + 1] << 16);
    }
    else
    {
        for (i = 0; i + 1 < an; i += 2)
            a2[i / 2] = a[i] | (a[i + 1] << 17);

        if (!squaring)
            for (i = 0; i + 1 < bn; i += 2)
                b2[i / 2] = b[i] | (b[i + 1] << 17);
    }

    if (an % 2) a2[an / 2] = a[an - 1];
    if (bn % 2) b2[bn / 2] = b[bn - 1];

    nmod_init(&mod2, FFT_PRIME);
    _nmod_poly_mul_mid_default_mpn_ctx(z2, 0, zn2, a2, an2, b2, bn2, mod2);

    if (M == M16)
    {
        UNPACK(M16);
    }
    else
    {
        UNPACK(M17);
    }

    TMP_END;

    return 1;
}

int
_nmod_poly_mullow_fft_small_repack(nn_ptr z, nn_srcptr a, slong an, nn_srcptr b, slong bn, slong zn, nmod_t mod)
{
    uint32_t aeven;
    uint32_t aodd;
    uint32_t beven;
    uint32_t bodd;
    ulong c0_bound;
    ulong c1_bound;
    ulong c2_bound;
    ulong n = mod.n;
    slong i;

    an = FLINT_MIN(an, zn);
    bn = FLINT_MIN(bn, zn);

    if (n > FFT_SMALL_MAX_REPACK || FLINT_MAX(an, bn) > (1 << 19))
        return 0;

    if (an + bn > 2.3 * repack_limit_tab[n])
        return 0;

    c0_bound = FLINT_MIN(an, bn) * (n - 1) * (n - 1);

    /* Coefficients are < 2^16, and (2^16)^3 < FFT_PRIME. */
    if (c0_bound < M16)
        return _nmod_poly_mullow_fft_small_repack_m(z, a, an, b, bn, zn, M16, mod);

    /*
    Since packing space is very limited, we try a Cauchy-Schwarz bound if
    the quick inf-norm bound fails: for random polynomials, a[k]*b[k]
    will be significantly smaller than (n-1)^2 on average when n is small.

    A = sum_{k=0}^{an-1} a[k] x^k
    B = sum_{k=0}^{bn-1} b[k] x^k

    a' = even coefficients
    a'' = odd coefficients  (possibly zero-padded to the same length as a')

    A2 = sum_{k=0}^{an2-1} (a'[k] + a''[k] M) x^k
    B2 = sum_{k=0}^{bn2-1} (b'[k] + b''[k] M) x^k

    A2 B2 = sum_{k=0}^{an2-1} sum_{l=0}^{bn2-1} (a'[k] + a''[k] M) (b'[l] + b''[l] M) x^{k+l}

    A2 B2 = sum_{k=0}^{an2-1} sum_{l=0}^{bn2-1} a'[k] b'[l] x^{k+l}  +
            sum_{k=0}^{an2-1} sum_{l=0}^{bn2-1} (a'[k] b''[l] + a''[k] b'[l]) M x^{k+l}
            sum_{k=0}^{an2-1} sum_{l=0}^{bn2-1} (a''[k] b''[l]) M^2 x^{k+l}

    Overflow cannot occur if the following inequalities are satisfied:

      ||a'||_2 ||b'||_2 < M
      ||a'||_2 ||b''||_2 + ||a''||_2 ||b'||_2 < M
      ||a''||_2 ||b''||_2 M^2 < FFT_PRIME

    Note that we can often choose M = 2^17 even though M^3 > FFT_PRIME,
    since the ||a''||_2 ||b''||_2 typically will be quite a bit smaller than M.

    TODO: these bounds are pessimistic if the product is very unbalanced;
    in that case, it would be better to split the larger polynomial into
    chunks.

    TODO: we could also combine a 1-norm bound with an inf-norm bound
    (Young's inequality). This would be faster to evaluate than
    Cauchy-Schwarz, but typically be less tight.
    */

    aeven = (an % 2) ? a[an - 1] : 0;
    aodd = 0;
    for (i = 0; i + 1 < an; i += 2)
    {
        aeven += (uint32_t) a[i] * (uint32_t) a[i];
        aodd  += (uint32_t) a[i + 1] * (uint32_t) a[i + 1];
    }

    aeven = sqrt(aeven) + 1;
    aodd = sqrt(aodd) + 1;

    if (a == b && an == bn)
    {
        beven = aeven;
        bodd = aodd;
    }
    else
    {
        beven = (bn % 2) ? b[bn - 1] : 0;
        bodd = 0;
        for (i = 0; i + 1 < bn; i += 2)
        {
            beven += (uint32_t) b[i] * (uint32_t) b[i];
            bodd  += (uint32_t) b[i + 1] * (uint32_t) b[i + 1];
        }

        beven = sqrt(beven) + 1;
        bodd = sqrt(bodd) + 1;
    }

    c0_bound = aeven * beven;
    c1_bound = aeven * bodd + aodd * beven;
    c2_bound = aodd * bodd;

    // flint_printf("n = %wu  an = %wd  bn = %wd  a' = %u  a'' = %u  b' = %u  b'' = %u   c0 %wu  c1 %wu  c2 %wu\n",
    //     n, an, bn, aeven, aodd, beven, bodd, c0_bound, c1_bound, c2_bound);

    if (c0_bound < M16 && c1_bound < M16 && c2_bound * M16 * M16 < FFT_PRIME)
        return _nmod_poly_mullow_fft_small_repack_m(z, a, an, b, bn, zn, M16, mod);

    if (c0_bound < M17 && c1_bound < M17 && c2_bound * M17 * M17 < FFT_PRIME)
        return _nmod_poly_mullow_fft_small_repack_m(z, a, an, b, bn, zn, M17, mod);

    return 0;
}

void
_nmod_poly_mullow_fft_small(nn_ptr z, nn_srcptr a, slong an, nn_srcptr b, slong bn, slong zn, nmod_t mod)
{
    if (!_nmod_poly_mullow_fft_small_repack(z, a, an, b, bn, zn, mod))
        _nmod_poly_mul_mid_default_mpn_ctx(z, 0, zn, a, an, b, bn, mod);
}


static const short fft_mul_tab[] = {1326, 1326, 1095, 802, 674, 537, 330, 306, 290,
274, 200, 192, 182, 173, 163, 99, 97, 93, 90, 82, 80, 438, 414, 324, 393,
298, 298, 268, 187, 185, 176, 176, 168, 167, 158, 158, 97, 96, 93, 92, 89,
89, 85, 85, 80, 81, 177, 172, 163, 162, 164, 176, 171, 167, 167, 164, 163,
163, 160, 165, 95, 96, 90, 94, };

static const short fft_sqr_tab[] = {1420, 1420, 1353, 964, 689, 569, 407, 353, 321,
321, 292, 279, 200, 182, 182, 159, 159, 152, 145, 139, 723, 626, 626, 569,
597, 448, 542, 292, 292, 200, 191, 191, 182, 182, 166, 166, 166, 159, 159,
159, 152, 152, 145, 145, 93, 200, 191, 182, 182, 182, 182, 191, 191, 191,
182, 182, 174, 182, 182, 182, 152, 152, 152, 145, };

/* todo: separate squaring table */
/* todo: check unbalanced cutoffs */
static const short fft_mullow_tab[] = {1115, 1115, 597, 569, 407, 321, 306, 279, 191,
182, 166, 159, 152, 145, 139, 89, 85, 78, 75, 75, 69, 174, 174, 166, 159,
152, 152, 152, 97, 101, 106, 111, 101, 101, 101, 139, 145, 145, 139, 145,
145, 139, 145, 145, 145, 182, 182, 182, 182, 182, 182, 191, 200, 220, 210,
200, 210, 210, 210, 210, 191, 182, 182, 174, };

int
_nmod_poly_mullow_want_fft_small(slong len1, slong len2, slong n, int squaring, nmod_t mod)
{
    slong len, bits, cutoff_len;

    if (len2 > len1)
        FLINT_SWAP(slong, len1, len2);

    if (n == len1 + len2 - 1)
    {
        if (mod.n <= FFT_SMALL_MAX_REPACK)
        {
            if (len2 < 64)
                return 0;

            len = len1 + len2 - 1;

            /* there is a big slowdown exactly at length 513
               as the transform length doubles */
            return len > 370 && !(len > 512 && len < 660)
                             && !(mod.n <= 3 && len > 1024 && len < 1100);
        }

        bits = NMOD_BITS(mod);
        cutoff_len = FLINT_MIN(len1, 2 * len2);

        if (squaring)
            return cutoff_len >= fft_sqr_tab[bits - 1];
        else
            return cutoff_len >= fft_mul_tab[bits - 1];
    }
    else
    {
        if (mod.n <= FFT_SMALL_MAX_REPACK)
        {
            if (len2 < 64)
                return 0;

            len = len1 + len2 - 1;
            return len > 450 && !(len > 512 && len < 700);
        }

        bits = NMOD_BITS(mod);
        cutoff_len = FLINT_MIN(len1, len2);

        return cutoff_len >= fft_mullow_tab[bits - 1];
    }
}

#else

int
_nmod_poly_mullow_fft_small_repack(nn_ptr z, nn_srcptr a, slong an, nn_srcptr b, slong bn, slong zn, nmod_t mod)
{
    return 0;
}

int
_nmod_poly_mullow_want_fft_small(slong len1, slong len2, slong n, int squaring, nmod_t mod)
{
    return 0;
}

void
_nmod_poly_mullow_fft_small(nn_ptr z, nn_srcptr a, slong an, nn_srcptr b, slong bn, slong zn, nmod_t mod)
{
    flint_throw(FLINT_ERROR, "fft_small is not available");
}


#endif

