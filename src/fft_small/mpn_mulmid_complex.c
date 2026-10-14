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
#include "fft_small.h"

/*
    Windowed complex products on transformed operands: the limbs
    [lo, hi) of the real and imaginary parts of (ar + i ai)(br + i bi),
    as two's-complement windows.

    The integers are split into chunks of P->bits bits and multiplied as
    chunk sequences, as in the real multiplication driver: four forward
    transforms, the bilinear combinations formed pointwise, two inverse
    transforms, and two signed exports of the windows. The windowing
    happens at the level of the chunk coefficients, also as in the real
    driver (mpn_mul.c): only the coefficients c_lo <= k < c_hi can reach
    the limb window, where c_lo is the first coefficient whose
    reconstruction (coeff_len limbs) reaches limb lo and c_hi the first
    lying entirely at or above limb hi. So

      - the exports reconstruct only those coefficients, and
      - when a power of two 2^d with max(operand chunks, c_hi) <= 2^d
        <= zlen and zlen - 2^d <= c_lo exists, the transforms are cyclic
        of length 2^d: the coefficients k >= 2^d alias onto k - 2^d < c_lo,
        which are never read as exact values.

    Error. The exports read whole blocks of slots, so some slots below
    c_lo -- dropped coefficients, or aliased sums of two -- are read.
    Every slot value is lifted to the centered range (-prod/2, prod/2]
    of the prime product, below 2^(FLINT_BITS coeff_len - 1) in absolute
    value, and a slot k < c_lo starts below bit FLINT_BITS (lo - coeff_len),
    so whatever it holds lies below limb lo; the slots below c_lo hold
    less than B^lo / 2 in total either way, both as computed and in the
    exact value. The window is therefore within one unit of
    floor(x / B^lo) for the exact part x, in either direction, and exact
    for lo = 0.

    The signs are applied in transform space, by choosing additions or
    subtractions and at most one pointwise negation per part, so that
    the exported windows are those of the actual parts (negating a
    floor window afterwards would cost another unit).
*/

static const ulong _zero_limb[1] = { 0 };

/* operand lengths: normalized, at least one (zero) limb */
static void
_part(nn_srcptr * x, slong * xn, slong hi)
{
    slong n = FLINT_MIN(*xn, hi);

    while (n > 0 && (*x)[n - 1] == 0)
        n--;

    if (n == 0)
    {
        *x = _zero_limb;
        n = 1;
    }

    *xn = n;
}

/* an mpn plan for signed two-term accumulations of products of operands
   of an_max x bn_max limbs, rewindowed to the chunk coefficients that
   can reach limbs [lo, hi), with the cyclic transform when it applies */
static int
_plan_window(fft_small_plan_t P, mpn_ctx_struct * R, slong an_max,
             slong bn_max, slong terms, slong lo, slong hi)
{
    ulong bits, nn, alen, blen, zlen, atrunc, btrunc, c_lo, c_hi, d, w;

    if (!fft_small_plan_init_mpn(P, R, an_max, bn_max, terms, 1))
        return 0;

    bits = P->bits;
    nn = (P->crts + P->np - 1)->coeff_len;
    alen = n_cdiv(FLINT_BITS * (ulong) an_max, bits);
    blen = n_cdiv(FLINT_BITS * (ulong) bn_max, bits);
    zlen = alen + blen - 1;
    atrunc = n_round_up(alen, BLK_SZ);
    btrunc = n_round_up(blen, BLK_SZ);

    c_lo = ((ulong) lo > nn) ? n_min(zlen, n_cdiv(((ulong) lo - nn) * FLINT_BITS, bits)) : 0;
    c_hi = n_min(zlen, n_cdiv((ulong) hi * FLINT_BITS, bits));

    d = n_flog2(zlen);
    w = n_pow2(d);

    if (d >= LG_BLK_SZ && d < P->depth &&
        atrunc <= w && btrunc <= w && c_hi <= w && w <= zlen && zlen <= c_lo + w)
    {
        P->depth = d;
        P->ztrunc = w;
        P->zl = c_lo;
        P->zh = c_hi;
        P->zn = zlen;

        /* the two-prime plans defer their normalization to the export,
           which indexes its scale factors by depth */
        if (P->np == 2)
            P->stride = sd_fft_ctx_data_size(d);
        else
            _fft_small_plan_set_normalizers(P);
    }

    return 1;
}

/* the number of slots the inverse transforms computed */
static ulong
_plan_nslots(const fft_small_plan_t P, slong an_max, slong bn_max)
{
    ulong alen = n_cdiv(FLINT_BITS * (ulong) an_max, P->bits);
    ulong blen = n_cdiv(FLINT_BITS * (ulong) bn_max, P->bits);
    return n_min(P->ztrunc, alen + blen - 1);
}

/* storage for k operands: a reservation in the context's persistent
   buffer, so that the transforms run on warm memory (fresh pages cost a
   tenth of the time at a few hundred kB per operand), or a plain
   allocation when another reservation is live */
static double *
_ops_storage(mpn_ctx_struct * R, const fft_small_plan_t P, ulong k, int * reserved)
{
    ulong sz = fft_small_op_sizeof_data(P);
    double * d = (double *) mpn_ctx_fit_buffer_reserve(R, 0, k * sz);

    *reserved = (d != NULL);
    if (d == NULL)
        d = flint_aligned_alloc(FLINT_FFT_SMALL_ALIGNMENT, k * sz);

    return d;
}

static void
_ops_storage_clear(mpn_ctx_struct * R, double * d, int reserved)
{
    if (reserved)
        mpn_ctx_fit_buffer_release(R);
    else
        flint_aligned_free(d);
}

/* Z = (-1)^e1 X1 Y1 + (-1)^e2 X2 Y2; Z may be X1 */
static void
_signed_pair(fft_small_op_t Z, const fft_small_op_t X1, const fft_small_op_t Y1,
             int e1, const fft_small_op_t X2, const fft_small_op_t Y2, int e2,
             const fft_small_plan_t P)
{
    fft_small_op_mul(Z, X1, Y1, P);
    if (e1 == e2)
        fft_small_op_addmul(Z, X2, Y2, P);
    else
        fft_small_op_submul(Z, X2, Y2, P);
    if (e1)
        fft_small_op_neg(Z, Z, P);
}

/*
    A zero part: every output part is then (at most) a single product of
    parts, re = ar br or -ai bi and im = ar bi or ai br, so only the
    nonzero parts are transformed -- an operand with one nonzero part
    once for both outputs -- and the plan holds one product per slot.
    The signs are still applied in transform space.
*/
static int
_mulmid_complex_sparse(nn_ptr zr, nn_ptr zi, slong lo, slong hi,
                       nn_srcptr * x, const slong * xn, const int * sgn,
                       const int * zero)
{
    /* output k: the product of parts _sa[k][j] and _sb[k][j], negated
       if _sn[k][j]; parts 0 ar, 1 ai, 2 br, 3 bi */
    static const int _sa[2][2] = { { 0, 1 }, { 0, 1 } };
    static const int _sb[2][2] = { { 2, 3 }, { 3, 2 } };
    static const int _sn[2][2] = { { 0, 1 }, { 0, 0 } };
    mpn_ctx_struct * R = get_default_mpn_ctx();
    fft_small_plan_t P;
    fft_small_op_struct X[6];
    nn_ptr z[2];
    slong an_max, bn_max;
    ulong nslots, sz, nops, i;
    int j, k, t[2];
    double * data;
    int reserved;

    z[0] = zr;
    z[1] = zi;

    for (k = 0; k < 2; k++)
    {
        t[k] = -1;
        for (j = 0; j < 2; j++)
            if (!zero[_sa[k][j]] && !zero[_sb[k][j]])
                t[k] = j;     /* at most one, as some part is zero */
    }

    if (t[0] < 0 && t[1] < 0)
    {
        flint_mpn_zero(zr, hi - lo);
        flint_mpn_zero(zi, hi - lo);
        return 1;
    }

    an_max = FLINT_MAX(xn[0], xn[1]);
    bn_max = FLINT_MAX(xn[2], xn[3]);

    if (!_plan_window(P, R, an_max, bn_max, 1, lo, hi))
        return 0;

    nslots = _plan_nslots(P, an_max, bn_max);

    nops = 2;
    for (j = 0; j < 4; j++)
        nops += !zero[j];

    data = _ops_storage(R, P, nops, &reserved);
    sz = fft_small_op_sizeof_data(P) / sizeof(double);

    i = 0;
    for (j = 0; j < 4; j++)
        if (!zero[j])
        {
            fft_small_op_init_borrowed(X + j, P, data + (i++) * sz);
            fft_small_fft_mpn(X + j, x[j], xn[j], P);
        }
    fft_small_op_init_borrowed(X + 4, P, data + (i++) * sz);
    fft_small_op_init_borrowed(X + 5, P, data + (i++) * sz);

    for (k = 0; k < 2; k++)
    {
        int a, b;

        if (t[k] < 0)
        {
            flint_mpn_zero(z[k], hi - lo);
            continue;
        }

        a = _sa[k][t[k]];
        b = _sb[k][t[k]];
        fft_small_op_mul(X + 4 + k, X + a, X + b, P);
        if (sgn[a] ^ sgn[b] ^ _sn[k][t[k]])
            fft_small_op_neg(X + 4 + k, X + 4 + k, P);
        fft_small_ifft(X + 4 + k, P);
        fft_small_export_mpn_signed_window(z[k], hi - lo, X + 4 + k, nslots, lo, P);
    }

    for (j = 0; j < 4; j++)
        if (!zero[j])
            fft_small_op_clear(X + j);
    fft_small_op_clear(X + 4);
    fft_small_op_clear(X + 5);
    _ops_storage_clear(R, data, reserved);
    fft_small_plan_clear(P);

    return 1;
}

int
_fft_small_mpn_mulmid_complex(nn_ptr zr, nn_ptr zi, slong lo, slong hi,
    nn_srcptr ar, slong arn, int ar_sgn,
    nn_srcptr ai, slong ain, int ai_sgn,
    nn_srcptr br, slong brn, int br_sgn,
    nn_srcptr bi, slong bin, int bi_sgn)
{
    mpn_ctx_struct * R = get_default_mpn_ctx();
    fft_small_plan_t P;
    fft_small_op_t Ar, Ai, Br, Bi, Z;
    slong an_max, bn_max;
    ulong nslots, sz;
    double * data;
    int reserved;

    FLINT_ASSERT(0 <= lo && lo < hi);

    _part(&ar, &arn, hi);
    _part(&ai, &ain, hi);
    _part(&br, &brn, hi);
    _part(&bi, &bin, hi);

    if (arn == 1 || ain == 1 || brn == 1 || bin == 1)
    {
        nn_srcptr x[4];
        slong xn[4];
        int sgn[4], zero[4], j;

        x[0] = ar; x[1] = ai; x[2] = br; x[3] = bi;
        xn[0] = arn; xn[1] = ain; xn[2] = brn; xn[3] = bin;
        sgn[0] = ar_sgn; sgn[1] = ai_sgn; sgn[2] = br_sgn; sgn[3] = bi_sgn;

        for (j = 0; j < 4; j++)
            zero[j] = (xn[j] == 1 && x[j][0] == 0);

        if (zero[0] || zero[1] || zero[2] || zero[3])
            return _mulmid_complex_sparse(zr, zi, lo, hi, x, xn, sgn, zero);
    }

    an_max = FLINT_MAX(arn, ain);
    bn_max = FLINT_MAX(brn, bin);

    if (!_plan_window(P, R, an_max, bn_max, 2, lo, hi))
        return 0;

    nslots = _plan_nslots(P, an_max, bn_max);

    data = _ops_storage(R, P, 5, &reserved);
    sz = fft_small_op_sizeof_data(P) / sizeof(double);
    fft_small_op_init_borrowed(Ar, P, data + 0 * sz);
    fft_small_op_init_borrowed(Ai, P, data + 1 * sz);
    fft_small_op_init_borrowed(Br, P, data + 2 * sz);
    fft_small_op_init_borrowed(Bi, P, data + 3 * sz);
    fft_small_op_init_borrowed(Z, P, data + 4 * sz);

    fft_small_fft_mpn(Ar, ar, arn, P);
    fft_small_fft_mpn(Ai, ai, ain, P);
    fft_small_fft_mpn(Br, br, brn, P);
    fft_small_fft_mpn(Bi, bi, bin, P);

    /* re = (-1)^(ar_sgn + br_sgn) ar br + (-1)^(ai_sgn + bi_sgn + 1) ai bi */
    _signed_pair(Z, Ar, Br, ar_sgn ^ br_sgn, Ai, Bi, ai_sgn ^ bi_sgn ^ 1, P);
    /* im = (-1)^(ar_sgn + bi_sgn) ar bi + (-1)^(ai_sgn + br_sgn) ai br,
       accumulated over ar, which dies at its first use */
    _signed_pair(Ar, Ar, Bi, ar_sgn ^ bi_sgn, Ai, Br, ai_sgn ^ br_sgn, P);

    fft_small_ifft(Z, P);
    fft_small_export_mpn_signed_window(zr, hi - lo, Z, nslots, lo, P);
    fft_small_ifft(Ar, P);
    fft_small_export_mpn_signed_window(zi, hi - lo, Ar, nslots, lo, P);

    fft_small_op_clear(Ar);
    fft_small_op_clear(Ai);
    fft_small_op_clear(Br);
    fft_small_op_clear(Bi);
    fft_small_op_clear(Z);
    _ops_storage_clear(R, data, reserved);
    fft_small_plan_clear(P);

    return 1;
}

int
_fft_small_mpn_sqrmid_complex(nn_ptr zr, nn_ptr zi, slong lo, slong hi,
    nn_srcptr ar, slong arn, int ar_sgn,
    nn_srcptr ai, slong ain, int ai_sgn)
{
    mpn_ctx_struct * R = get_default_mpn_ctx();
    fft_small_plan_t P;
    fft_small_op_t Ar, Ai, Z;
    slong an_max;
    ulong nslots, sz;
    double * data;
    int reserved;

    FLINT_ASSERT(0 <= lo && lo < hi);

    _part(&ar, &arn, hi);
    _part(&ai, &ain, hi);
    an_max = FLINT_MAX(arn, ain);

    if (!_plan_window(P, R, an_max, an_max, 2, lo, hi))
        return 0;

    nslots = _plan_nslots(P, an_max, an_max);

    data = _ops_storage(R, P, 3, &reserved);
    sz = fft_small_op_sizeof_data(P) / sizeof(double);
    fft_small_op_init_borrowed(Ar, P, data + 0 * sz);
    fft_small_op_init_borrowed(Ai, P, data + 1 * sz);
    fft_small_op_init_borrowed(Z, P, data + 2 * sz);

    fft_small_fft_mpn(Ar, ar, arn, P);
    fft_small_fft_mpn(Ai, ai, ain, P);

    /* re = ar^2 - ai^2, sign free */
    fft_small_op_sqr(Z, Ar, P);
    fft_small_op_subsqr(Z, Ai, P);

    /* im = (-1)^(ar_sgn + ai_sgn) 2 ar ai */
    fft_small_op_mul(Ar, Ar, Ai, P);
    fft_small_op_add(Ar, Ar, Ar, P);
    if (ar_sgn ^ ai_sgn)
        fft_small_op_neg(Ar, Ar, P);

    fft_small_ifft(Z, P);
    fft_small_export_mpn_signed_window(zr, hi - lo, Z, nslots, lo, P);
    fft_small_ifft(Ar, P);
    fft_small_export_mpn_signed_window(zi, hi - lo, Ar, nslots, lo, P);

    fft_small_op_clear(Ar);
    fft_small_op_clear(Ai);
    fft_small_op_clear(Z);
    _ops_storage_clear(R, data, reserved);
    fft_small_plan_clear(P);

    return 1;
}
