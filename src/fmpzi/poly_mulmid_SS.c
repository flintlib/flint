/*
    Copyright (C) 2026 Fredrik Johansson

    Developed using Claude Opus 5.5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpzi.h"
#include "fft.h"
#include "fft_small.h"

/*
    Middle products of Gaussian integer polynomials by Schoenhage-Strassen
    convolution.

    In the coefficient ring Z/(2^N + 1) of the transform, s = 2^(N/2)
    satisfies s^2 = -1. Since 2^N + 1 = (s + i)(s - i) with coprime
    factors in Z[i],

        Z[i]/(2^N + 1)  =  Z/(2^N + 1) x Z/(2^N + 1),
        a + bi         ->  (a + s b, a - s b),

    and both components are ring homomorphisms applied coefficientwise.
    A Gaussian convolution is therefore two ordinary convolutions of the
    same length and ring size: of u = re + s im and of v = re - s im.
    The results give back

        re = (U + V) / 2,    im = (V - U) s / 2      (mod 2^N + 1),

    exactly once N exceeds the bit size of the components of the result.
    All the conversions are shifts and additions. The cost is six
    transforms and two pointwise products per point, against six
    transforms and three or four pointwise products for real-arithmetic
    schemes sharing transforms, and nine or twelve transforms for three
    or four independent real products.

    The convolution length is chosen as for _fmpz_poly_mulmid_SS,
    including the cyclic (wraparound) option for middle products.

    Two variants follow below: the single embedding a + bi -> a + s b in
    one convolution over a ring of twice the size, and the products of
    a real or imaginary operand, which need no embedding.
*/

/* u <- u + s v, v <- u - s v (u holds the real part, v the imaginary) */
static void
_embed_pair(nn_ptr u, nn_ptr v, slong limbs, nn_ptr tmp)
{
    fft_adjust(tmp, v, FLINT_BITS * limbs / 2, limbs, 1);
    butterfly_lshB(u, v, u, tmp, limbs, 0, 0);
    mpn_normmod_2expp1(u, limbs);
    mpn_normmod_2expp1(v, limbs);
}

/* re <- (U + V)/2, im <- (V - U) s/2 */
static void
_unembed_pair(nn_ptr re, nn_ptr im, nn_srcptr U, nn_srcptr V,
              slong limbs, nn_ptr tmp1, nn_ptr tmp2)
{
    butterfly_lshB(tmp1, tmp2, (nn_ptr) V, (nn_ptr) U, limbs, 0, 0);
    mpn_div_2expmod_2expp1(re, tmp1, limbs, 1);
    mpn_normmod_2expp1(re, limbs);
    fft_adjust(tmp1, tmp2, FLINT_BITS * limbs / 2, limbs, 1);
    mpn_div_2expmod_2expp1(im, tmp1, limbs, 1);
    mpn_normmod_2expp1(im, limbs);
}

static slong
_fmpzi_vec_max_bits_abs(const fmpzi_struct * x, slong len)
{
    slong i, b = 0;

    for (i = 0; i < len; i++)
    {
        slong br = fmpz_bits(fmpzi_realref(x + i));
        slong bi = fmpz_bits(fmpzi_imagref(x + i));
        b = FLINT_MAX(b, br);
        b = FLINT_MAX(b, bi);
    }

    return b;
}

/* pointer array of length 4n into a fresh block of 4n * size limbs */
static nn_ptr *
_fft_vec_alloc(slong n, slong size)
{
    slong i;
    nn_ptr * v, ptr;

    v = flint_malloc(4 * (n + n * size) * sizeof(ulong));
    for (i = 0, ptr = (nn_ptr) v + 4 * n; i < 4 * n; i++, ptr += size)
        v[i] = ptr;

    return v;
}

/* load the real and imaginary parts of (poly, len), embed, zero pad */
static void
_fft_vec_load_pair(nn_ptr * U, nn_ptr * V, const fmpzi_struct * poly,
                   slong len, slong n, slong limbs, nn_ptr tmp)
{
    slong i;

    for (i = 0; i < len; i++)
    {
        _fmpz_vec_get_fft(U + i, fmpzi_realref(poly + i), limbs, 1);
        _fmpz_vec_get_fft(V + i, fmpzi_imagref(poly + i), limbs, 1);
        _embed_pair(U[i], V[i], limbs, tmp);
    }

    for (i = len; i < 4 * n; i++)
    {
        flint_mpn_zero(U[i], limbs + 1);
        flint_mpn_zero(V[i], limbs + 1);
    }
}

/* 1 if every imaginary part is zero, 2 if every real part is (and some
   imaginary part not), 0 otherwise */
static int
_fmpzi_vec_single_part(const fmpzi_struct * x, slong len)
{
    slong i;
    int re = 1, im = 1;

    for (i = 0; i < len && (re || im); i++)
    {
        re = re && fmpz_is_zero(fmpzi_imagref(x + i));
        im = im && fmpz_is_zero(fmpzi_realref(x + i));
    }

    return re ? 1 : im ? 2 : 0;
}

/* load part k (1 real, 2 imaginary) of (poly, len), zero pad */
static void
_fft_vec_load_part(nn_ptr * U, const fmpzi_struct * poly, slong len, int k,
                   slong n, slong limbs)
{
    slong i;

    for (i = 0; i < len; i++)
        _fmpz_vec_get_fft(U + i, (k == 1) ? fmpzi_realref(poly + i)
                                          : fmpzi_imagref(poly + i), limbs, 1);

    for (i = len; i < 4 * n; i++)
        flint_mpn_zero(U[i], limbs + 1);
}

/*
    poly1 has a single nonzero part p (k1 = 1 real, 2 imaginary): the
    product is p times the parts of poly2, two real convolutions sharing
    the transform of p (precached), or one when poly2 also has a single
    part. Without sums of products the coefficients need one bit less
    than the pair embedding, which is not used.
*/
static void
_fmpzi_poly_mulmid_SS_part(fmpzi_struct * res,
    const fmpzi_struct * poly1, slong len1, int k1,
    const fmpzi_struct * poly2, slong len2, int k2,
    slong nlo, slong nhi, slong bits1, slong bits2, int squaring)
{
    slong len_out, conv_len, depth, n, limbs, size, negm, res_bits, i, N;
    nn_ptr * jj, * i1, * i2;
    nn_ptr * t1, * t2, * s1, * tt;
    nn_ptr ptr;
#if FLINT_HAVE_FFT_SMALL
    sd_fft_mpn_mulmod_2expp1_ctx_t sdctx;
#endif

    FLINT_ASSERT(k1 != 0);

    len_out = len1 + len2 - 1;
    depth = fft_convolution_choose_mid(&conv_len, len_out, nlo, nhi);
    n = WORD(1) << depth;

    /* |c| < min(len1, len2) 2^(bits1 + bits2), plus a sign bit */
    res_bits = bits1 + bits2 + FLINT_CLOG2(FLINT_MIN(len1, len2)) + 1;
    res_bits = (((res_bits - 1) >> depth) + 1) << depth;

    limbs = (res_bits - 1) / FLINT_BITS + 1;
    limbs = fft_adjust_limbs_sd_fft(limbs, n, &negm);
    size = limbs + 1;

    N = flint_get_num_threads();
    jj = _fft_vec_alloc(n, size);
    i1 = squaring ? jj : _fft_vec_alloc(n, size);
    i2 = (k2 == 0) ? _fft_vec_alloc(n, size) : NULL;

    ptr = flint_malloc(5 * size * N * sizeof(ulong));
    t1 = flint_malloc(4 * N * sizeof(nn_ptr));
    t2 = t1 + N;
    s1 = t2 + N;
    tt = s1 + N;

    t1[0] = ptr;
    t2[0] = t1[0] + size * N;
    s1[0] = t2[0] + size * N;
    tt[0] = s1[0] + size * N;
    for (i = 1; i < N; i++)
    {
        t1[i] = t1[i - 1] + size;
        t2[i] = t2[i - 1] + size;
        s1[i] = s1[i - 1] + size;
        tt[i] = tt[i - 1] + 2 * size;
    }

    _fft_vec_load_part(jj, poly1, len1, k1, n, limbs);
    if (k2 == 0)
    {
        _fft_vec_load_part(i1, poly2, len2, 1, n, limbs);
        _fft_vec_load_part(i2, poly2, len2, 2, n, limbs);
    }
    else if (!squaring)
        _fft_vec_load_part(i1, poly2, len2, k2, n, limbs);

#if FLINT_HAVE_FFT_SMALL
    if (negm != 0)
        sd_fft_mpn_mulmod_2expp1_ctx_init(sdctx, get_default_mpn_ctx(),
                                          FLINT_BITS * limbs, negm);
#endif

    if (k2 == 0)
    {
        /* one forward transform of p for both products */
        fft_precache(jj, depth, limbs, conv_len, t1, t2, s1);
#if FLINT_HAVE_FFT_SMALL
        if (negm != 0)
        {
            fft_convolution_precache_sd_fft(sdctx, i1, jj, depth, limbs,
                                            conv_len, t1, t2, s1, tt);
            fft_convolution_precache_sd_fft(sdctx, i2, jj, depth, limbs,
                                            conv_len, t1, t2, s1, tt);
        }
        else
#endif
        {
            fft_convolution_precache(i1, jj, depth, limbs, conv_len, t1, t2, s1, tt);
            fft_convolution_precache(i2, jj, depth, limbs, conv_len, t1, t2, s1, tt);
        }
    }
    else
    {
#if FLINT_HAVE_FFT_SMALL
        if (negm != 0)
            fft_convolution_sd_fft(sdctx, i1, jj, depth, limbs, conv_len,
                                   t1, t2, s1, tt);
        else
#endif
            fft_convolution(i1, jj, depth, limbs, conv_len, t1, t2, s1, tt);
    }

#if FLINT_HAVE_FFT_SMALL
    if (negm != 0)
        sd_fft_mpn_mulmod_2expp1_ctx_clear(sdctx);
#endif

    /* res = u (C1 + i C2) for poly2 with both parts (C1 = p re, C2 =
       p im), res = u v C1 for a single part, with the units u, v = 1
       or i of the parts */
    for (i = nlo; i < nhi; i++)
    {
        fmpzi_struct * r = res + i - nlo;

        if (k2 == 0)
        {
            fmpz * c1 = (k1 == 1) ? fmpzi_realref(r) : fmpzi_imagref(r);
            fmpz * c2 = (k1 == 1) ? fmpzi_imagref(r) : fmpzi_realref(r);

            _fmpz_vec_set_fft(c1, 1, i1 + i, limbs, 1);
            _fmpz_vec_set_fft(c2, 1, i2 + i, limbs, 1);
            if (k1 == 2)
                fmpz_neg(c2, c2);
        }
        else
        {
            int k = k1 + k2;     /* 2: 1, 3: i, 4: -1 */
            fmpz * c = (k == 3) ? fmpzi_imagref(r) : fmpzi_realref(r);

            _fmpz_vec_set_fft(c, 1, i1 + i, limbs, 1);
            fmpz_zero((k == 3) ? fmpzi_realref(r) : fmpzi_imagref(r));
            if (k == 4)
                fmpz_neg(c, c);
        }
    }

    flint_free(jj);
    if (!squaring)
        flint_free(i1);
    if (k2 == 0)
        flint_free(i2);
    flint_free(ptr);
    flint_free(t1);
}

/*
    The common preparation: truncation of the inputs to the window,
    len1 >= len2, the bit sizes, and the products that need no
    Gaussian convolution (a zero operand, or a real or imaginary one).
    Returns 1 if the window has been computed.
*/
static int
_fmpzi_poly_mulmid_SS_prepare(fmpzi_struct * res,
    const fmpzi_struct ** p1, slong * l1, const fmpzi_struct ** p2,
    slong * l2, slong * lo, slong * hi, slong * b1, slong * b2, int * sq)
{
    const fmpzi_struct * poly1 = *p1, * poly2 = *p2;
    slong len1 = *l1, len2 = *l2, nlo = *lo, nhi = *hi;
    slong bits1, bits2, i;
    int squaring;

    FLINT_ASSERT(len1 != 0);
    FLINT_ASSERT(len2 != 0);
    FLINT_ASSERT(nlo >= 0);
    FLINT_ASSERT(nlo < nhi);
    FLINT_ASSERT(nhi <= len1 + len2 - 1);

    /* low truncation of the inputs */
    len1 = FLINT_MIN(len1, nhi);
    len2 = FLINT_MIN(len2, nhi);

    /* high truncation of the inputs */
    {
        slong nlo2 = (len1 + len2 - 1) - nlo;

        if (len1 > nlo2)
        {
            slong trunc = len1 - nlo2;
            poly1 += trunc;
            len1 -= trunc;
            nlo -= trunc;
            nhi -= trunc;
        }

        if (len2 > nlo2)
        {
            slong trunc = len2 - nlo2;
            poly2 += trunc;
            len2 -= trunc;
            nlo -= trunc;
            nhi -= trunc;
        }
    }

    if (len1 < len2)
    {
        FLINT_SWAP(const fmpzi_struct *, poly1, poly2);
        FLINT_SWAP(slong, len1, len2);
    }

    squaring = (poly1 == poly2 && len1 == len2);

    bits1 = _fmpzi_vec_max_bits_abs(poly1, len1);
    bits2 = squaring ? bits1 : _fmpzi_vec_max_bits_abs(poly2, len2);

    if (bits1 == 0 || bits2 == 0)
    {
        for (i = 0; i < nhi - nlo; i++)
            fmpzi_zero(res + i);
        return 1;
    }

    /* a real or imaginary operand */
    {
        int k1 = _fmpzi_vec_single_part(poly1, len1);
        int k2 = squaring ? k1 : _fmpzi_vec_single_part(poly2, len2);

        if (k1 != 0 || k2 != 0)
        {
            if (k1 == 0)
                _fmpzi_poly_mulmid_SS_part(res, poly2, len2, k2, poly1,
                                     len1, k1, nlo, nhi, bits2, bits1, 0);
            else
                _fmpzi_poly_mulmid_SS_part(res, poly1, len1, k1, poly2,
                                     len2, k2, nlo, nhi, bits1, bits2, squaring);
            return 1;
        }
    }

    *p1 = poly1; *l1 = len1; *p2 = poly2; *l2 = len2;
    *lo = nlo; *hi = nhi; *b1 = bits1; *b2 = bits2; *sq = squaring;
    return 0;
}

/* bits of |re|, |im| < 2 len2 2^(bits1 + bits2) in the result, plus
   a sign bit (len2 the shorter length) */
static slong
_result_bits(slong bits1, slong bits2, slong len2)
{
    return bits1 + bits2 + FLINT_CLOG2(len2) + 1 + 1;
}

/* the ring of the pair method: the result bits, rounded up for the
   sqrt2 trick */
static slong
_pair_limbs(slong * negm, slong bits1, slong bits2, slong len2, slong depth)
{
    slong res_bits = _result_bits(bits1, bits2, len2), limbs;

    res_bits = (((res_bits - 1) >> depth) + 1) << depth;
    limbs = (res_bits - 1) / FLINT_BITS + 1;

    return fft_adjust_limbs_sd_fft(limbs, WORD(1) << depth, negm);
}

/* the ring of the single embedding: an even number of limbs, half of
   which hold one more bit than a result part (see below) */
static slong
_single_limbs(slong * negm, slong bits1, slong bits2, slong len2, slong depth)
{
    slong unit = FLINT_MAX(2 * FLINT_BITS, WORD(1) << depth);
    slong half = (_result_bits(bits1, bits2, len2) + 1 + FLINT_BITS - 1) / FLINT_BITS;
    slong bits = 2 * FLINT_BITS * half, limbs;

    for (;;)
    {
        bits = ((bits + unit - 1) / unit) * unit;
        limbs = fft_adjust_limbs_sd_fft(bits / FLINT_BITS, WORD(1) << depth, negm);
        if (limbs % 2 == 0)
            return limbs;
        bits = FLINT_BITS * (limbs + 1);
    }
}

static void
_SS_pair(fmpzi_struct * res,
    const fmpzi_struct * poly1, slong len1,
    const fmpzi_struct * poly2, slong len2,
    slong nlo, slong nhi, slong bits1, slong bits2, int squaring)
{
    slong len_out, conv_len, depth, n, limbs, size, negm;
    slong i, N;
    nn_ptr * u1, * v1, * u2, * v2;
    nn_ptr * t1, * t2, * s1, * tt;
    nn_ptr ptr, tmp;

    len_out = len1 + len2 - 1;
    depth = fft_convolution_choose_mid(&conv_len, len_out, nlo, nhi);
    n = WORD(1) << depth;

    limbs = _pair_limbs(&negm, bits1, bits2, len2, depth);
    size = limbs + 1;

    /* the transforms exchange pointers between the vectors and the
       temporaries, so everything stays allocated until the outputs
       have been read */
    N = flint_get_num_threads();
    u1 = _fft_vec_alloc(n, size);
    v1 = _fft_vec_alloc(n, size);
    u2 = squaring ? u1 : _fft_vec_alloc(n, size);
    v2 = squaring ? v1 : _fft_vec_alloc(n, size);

    ptr = flint_malloc((5 * size * N + 2 * size) * sizeof(ulong));
    t1 = flint_malloc(4 * N * sizeof(nn_ptr));
    t2 = t1 + N;
    s1 = t2 + N;
    tt = s1 + N;

    t1[0] = ptr;
    t2[0] = t1[0] + size * N;
    s1[0] = t2[0] + size * N;
    tt[0] = s1[0] + size * N;
    for (i = 1; i < N; i++)
    {
        t1[i] = t1[i - 1] + size;
        t2[i] = t2[i - 1] + size;
        s1[i] = s1[i - 1] + size;
        tt[i] = tt[i - 1] + 2 * size;
    }
    tmp = tt[0] + 2 * size * N;

    _fft_vec_load_pair(u1, v1, poly1, len1, n, limbs, tmp);
    if (!squaring)
        _fft_vec_load_pair(u2, v2, poly2, len2, n, limbs, tmp);

#if FLINT_HAVE_FFT_SMALL
    if (negm != 0)
    {
        sd_fft_mpn_mulmod_2expp1_ctx_t sdctx;
        sd_fft_mpn_mulmod_2expp1_ctx_init(sdctx, get_default_mpn_ctx(),
                                          FLINT_BITS * limbs, negm);
        fft_convolution_sd_fft(sdctx, u1, u2, depth, limbs, conv_len,
                               t1, t2, s1, tt);
        fft_convolution_sd_fft(sdctx, v1, v2, depth, limbs, conv_len,
                               t1, t2, s1, tt);
        sd_fft_mpn_mulmod_2expp1_ctx_clear(sdctx);
    }
    else
#endif
    {
        fft_convolution(u1, u2, depth, limbs, conv_len, t1, t2, s1, tt);
        fft_convolution(v1, v2, depth, limbs, conv_len, t1, t2, s1, tt);
    }

    /* recover the window. The convolutions only permute buffers among
       their vectors and the temporaries, so the buffers now held by t1,
       t2 and s1 are disjoint from those of u1 and v1 and serve as
       scratch. The inputs have been consumed, so res may alias them. */
    {
        nn_ptr re = t1[0], im = t2[0];

        for (i = nlo; i < nhi; i++)
        {
            _unembed_pair(re, im, u1[i], v1[i], limbs, s1[0], tmp);
            _fmpz_vec_set_fft(fmpzi_realref(res + i - nlo), 1, &re, limbs, 1);
            _fmpz_vec_set_fft(fmpzi_imagref(res + i - nlo), 1, &im, limbs, 1);
        }
    }

    flint_free(u1);
    flint_free(v1);
    if (!squaring)
    {
        flint_free(u2);
        flint_free(v2);
    }
    flint_free(ptr);
    flint_free(t1);
}

/*
    The single embedding: i -> s = 2^(N/2) maps a Gaussian polynomial
    to a + s b over Z/(2^N + 1) as a ring homomorphism, so one
    convolution gives T = re + s im (mod 2^N + 1) for the product. With
    N = 2 FLINT_BITS h and |re|, |im| < 2^(FLINT_BITS h - 2), the
    integer T lies in (-2^(N - 1), 2^(N - 1)), so both parts are read
    off its two's complement representation: re is the signed low half
    and im the signed high half plus the borrow (the top bit of the low
    half). One convolution at twice the ring size against two for the
    pair method: the same transform work and memory, but rings are
    rounded up to multiples of the transform length (and to sizes the
    pointwise products support) once instead of twice, so it wins when
    the coefficients are small for the length.
*/

/* load a + s b, as an (N + FLINT_BITS)-bit two's complement integer */
static void
_fft_vec_load_single(nn_ptr * U, const fmpzi_struct * poly, slong len,
                     slong n, slong limbs, nn_ptr tmp)
{
    slong i, h = limbs / 2;

    for (i = 0; i < len; i++)
    {
        fmpz_get_signed_ui_array(U[i], limbs + 1, fmpzi_realref(poly + i));
        fmpz_get_signed_ui_array(tmp, limbs + 1 - h, fmpzi_imagref(poly + i));
        mpn_add_n(U[i] + h, U[i] + h, tmp, limbs + 1 - h);
    }

    for (i = len; i < 4 * n; i++)
        flint_mpn_zero(U[i], limbs + 1);
}

/* the parts of a normalized residue X in [0, 2^N]; X is destroyed */
static void
_unembed_single(fmpzi_t res, nn_ptr X, slong limbs)
{
    slong h = limbs / 2;

    /* T < 0 iff X > 2^(N - 1); then X - 1 is the N-bit two's
       complement of T */
    if (X[limbs] != 0 || (X[limbs - 1] >> (FLINT_BITS - 1)) != 0)
        mpn_sub_1(X, X, limbs + 1, 1);

    fmpz_set_signed_ui_array(fmpzi_realref(res), X, h);
    fmpz_set_signed_ui_array(fmpzi_imagref(res), X + h, limbs - h);
    if ((X[h - 1] >> (FLINT_BITS - 1)) != 0)
        fmpz_add_ui(fmpzi_imagref(res), fmpzi_imagref(res), 1);
}

static void
_SS_single(fmpzi_struct * res,
    const fmpzi_struct * poly1, slong len1,
    const fmpzi_struct * poly2, slong len2,
    slong nlo, slong nhi, slong bits1, slong bits2, int squaring)
{
    slong len_out, conv_len, depth, n, limbs, size, negm;
    slong i, N;
    nn_ptr * u1, * u2;
    nn_ptr * t1, * t2, * s1, * tt;
    nn_ptr ptr, tmp;

    len_out = len1 + len2 - 1;
    depth = fft_convolution_choose_mid(&conv_len, len_out, nlo, nhi);
    n = WORD(1) << depth;

    limbs = _single_limbs(&negm, bits1, bits2, len2, depth);
    size = limbs + 1;

    FLINT_ASSERT((FLINT_BITS * limbs) % n == 0);

    N = flint_get_num_threads();
    u1 = _fft_vec_alloc(n, size);
    u2 = squaring ? u1 : _fft_vec_alloc(n, size);

    ptr = flint_malloc((5 * size * N + size) * sizeof(ulong));
    t1 = flint_malloc(4 * N * sizeof(nn_ptr));
    t2 = t1 + N;
    s1 = t2 + N;
    tt = s1 + N;

    t1[0] = ptr;
    t2[0] = t1[0] + size * N;
    s1[0] = t2[0] + size * N;
    tt[0] = s1[0] + size * N;
    for (i = 1; i < N; i++)
    {
        t1[i] = t1[i - 1] + size;
        t2[i] = t2[i - 1] + size;
        s1[i] = s1[i - 1] + size;
        tt[i] = tt[i - 1] + 2 * size;
    }
    tmp = tt[0] + 2 * size * N;

    _fft_vec_load_single(u1, poly1, len1, n, limbs, tmp);
    if (!squaring)
        _fft_vec_load_single(u2, poly2, len2, n, limbs, tmp);

#if FLINT_HAVE_FFT_SMALL
    if (negm != 0)
    {
        sd_fft_mpn_mulmod_2expp1_ctx_t sdctx;
        sd_fft_mpn_mulmod_2expp1_ctx_init(sdctx, get_default_mpn_ctx(),
                                          FLINT_BITS * limbs, negm);
        fft_convolution_sd_fft(sdctx, u1, u2, depth, limbs, conv_len,
                               t1, t2, s1, tt);
        sd_fft_mpn_mulmod_2expp1_ctx_clear(sdctx);
    }
    else
#endif
        fft_convolution(u1, u2, depth, limbs, conv_len, t1, t2, s1, tt);

    /* the inputs have been consumed, so res may alias them */
    for (i = nlo; i < nhi; i++)
        _unembed_single(res + i - nlo, u1[i], limbs);

    flint_free(u1);
    if (!squaring)
        flint_free(u2);
    flint_free(ptr);
    flint_free(t1);
}

void
_fmpzi_poly_mulmid_SS_pair(fmpzi_struct * res,
    const fmpzi_struct * poly1, slong len1,
    const fmpzi_struct * poly2, slong len2,
    slong nlo, slong nhi)
{
    slong bits1, bits2;
    int squaring;

    if (!_fmpzi_poly_mulmid_SS_prepare(res, &poly1, &len1, &poly2, &len2,
                                       &nlo, &nhi, &bits1, &bits2, &squaring))
        _SS_pair(res, poly1, len1, poly2, len2, nlo, nhi, bits1, bits2, squaring);
}

void
_fmpzi_poly_mulmid_SS_single(fmpzi_struct * res,
    const fmpzi_struct * poly1, slong len1,
    const fmpzi_struct * poly2, slong len2,
    slong nlo, slong nhi)
{
    slong bits1, bits2;
    int squaring;

    if (!_fmpzi_poly_mulmid_SS_prepare(res, &poly1, &len1, &poly2, &len2,
                                       &nlo, &nhi, &bits1, &bits2, &squaring))
        _SS_single(res, poly1, len1, poly2, len2, nlo, nhi, bits1, bits2, squaring);
}

/*
    Choice of method. The pair method runs six transforms of rings of
    lp limbs and two pointwise products, the single embedding three
    transforms of rings of ls limbs and one pointwise product, where
    ls is about 2 lp before both are rounded to the transform length
    and to the sizes the pointwise products support. Measured (x86-64):

      - pair rings up to 1280 bits: the single embedding, 15-35% faster
        (half the pointwise products and conversions of small rings);
      - pair rings from 32768 bits, where both pointwise products run on
        the negacyclic fft_small engine and cost about linearly: equal
        within a few percent, slightly favoring the single embedding
        when ls <= 2 lp;
      - in between, the pointwise products of the twice larger rings
        cost more than twice as much (up to twice as much in total at
        2048-bit parts), and the single embedding only wins when its
        ring is rounded at least 7% less, ls <= 0.93 (2 lp), which
        happens for small coefficients relative to the length.
*/
static int
_SS_use_single(slong lp, slong ls)
{
    if (FLINT_BITS * lp <= 1280)
        return 1;
    if (FLINT_BITS * lp >= 32768)
        return ls <= 2 * lp;
    return 100 * ls <= 186 * lp;
}

void
_fmpzi_poly_mulmid_SS(fmpzi_struct * res,
    const fmpzi_struct * poly1, slong len1,
    const fmpzi_struct * poly2, slong len2,
    slong nlo, slong nhi)
{
    slong bits1, bits2, conv_len, depth, lp, ls, negm;
    int squaring;

    if (_fmpzi_poly_mulmid_SS_prepare(res, &poly1, &len1, &poly2, &len2,
                                      &nlo, &nhi, &bits1, &bits2, &squaring))
        return;

    depth = fft_convolution_choose_mid(&conv_len, len1 + len2 - 1, nlo, nhi);
    lp = _pair_limbs(&negm, bits1, bits2, len2, depth);
    ls = _single_limbs(&negm, bits1, bits2, len2, depth);

    if (_SS_use_single(lp, ls))
        _SS_single(res, poly1, len1, poly2, len2, nlo, nhi, bits1, bits2, squaring);
    else
        _SS_pair(res, poly1, len1, poly2, len2, nlo, nhi, bits1, bits2, squaring);
}
