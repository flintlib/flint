/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "gmpcompat.h"
#include "fmpz.h"
#include "mpn_extras.h"
#include "ulong_extras.h"
#include "gr.h"

#if FLINT_HAVE_FFT_SMALL

#include "machine_vectors.h"
#include "fft_small.h"

/*
    Transformed big integers (fft_small), for bilinear expressions such as
    complex multiplications and matrix products over Z.

    A context has a capacity: results of expressions must have absolute
    value below 2^bits_bound, accumulate at most terms_bound elementary
    products per bit-chunk, and have multiplicative depth at most 2 in
    converted-in operands. Elements store the fft_small transform of the
    operand's bit-chunk sequence (the existing packing code, untouched)
    together with a virtual chunk length, an accumulation counter, a depth
    counter, a sign bit and a chunks-possibly-signed flag: an element
    represents (-1)^sign * (sum of chunks * 2^(j*bits)).

    Signs: conversion in takes a magnitude and a sign bit, so the packing
    never sees signed data. Multiplication multiplies sign bits and keeps
    chunks nonnegative. Additive operations compare sign bits and switch
    the roles of the pointwise additions and subtractions (add <-> sub,
    addmul <-> submul); when a subtraction actually occurs, the chunk
    values may become signed and the element is flagged. Conversion out of
    a flagged element biases every chunk by C = terms_bound * (2^bits-1)^2
    in the pointwise domain -- so the existing unsigned reconstruction and
    carry propagation run untouched -- which adds the known integer
    C * (1 + 2^bits + ... + 2^((m-1) bits)) to the result; one comparison
    against that bias integer determines the result's sign and one
    subtraction recovers its magnitude. Unsigned contexts and provably
    nonnegative elements skip the bias entirely.

    A context may be shared by any number of threads: after
    construction it only changes through its operand pool, whose
    bookkeeping is locked, and conversions take independent scratch. A
    single element belongs to one thread at a time. A context holding a
    reservation (the fit_buffer strategy) pins the creating thread's
    fft_small buffer and must be cleared on that thread.

    Windowed contexts (gr_ctx_init_transformed_mpn_window) plan for
    products of which only the limbs [lo, hi) are wanted, with cyclic
    transforms when the window allows; their products whose slots alias
    then convert only through gr_transformed_mpn_get_window inside that
    window, and the other conversions decline.
*/

typedef struct
{
    fft_small_plan_t P;
    slong bits_bound;               /* result magnitude < 2^bits_bound */
    slong zcap;                     /* usable product chunks (plan chunk
                                       range minus bias-digit headroom) */
    ulong terms_bound;
    ulong max_depth;
    int is_signed;
    ulong m_orig[MPN_CTX_NCRTS];    /* saved export normalizers */
    /* operand storage (the shared fft_small pool: a reservation in the
       creating thread's mpn_ctx buffer for the fit_buffer strategy, or
       nothing, elements then allocating their own data) and conversion
       staging (the pool's staging region or a dead operand slab, else a
       plain allocation per conversion). The pool is thread safe, so
       contexts may be shared between threads. */
    int alloc_strategy;
    int scratch_from_slab;
    fft_small_op_pool_t pool;
    /* the output window of a windowed context (limbs [win_lo, win_hi)),
       which a cyclic plan restricts the conversions of products to */
    int windowed;
    slong win_lo;
    slong win_hi;
} tmpn_ctx_struct;

typedef struct
{
    fft_small_op_struct op;
    slong nchunks;                  /* virtual chunk length; 0 = zero */
    ulong terms;
    ulong depth;
    int sign;                       /* 0 or 1: (-1)^sign */
    int negs;                       /* chunk values may be signed */
    slong small;                    /* small side value: the element
                                       represents (-1)^sign * T + small,
                                       with T = 0 when nchunks == 0.
                                       Coefficients 0 and +-1 are held
                                       here without any transform, and
                                       products of such coefficients
                                       accumulate here as integer
                                       additions, folded into the
                                       reconstruction at conversion
                                       out. */
} tmpn_struct;

#define TMPN_CTX(ctx) ((tmpn_ctx_struct *) GR_CTX_DATA_AS_PTR(ctx))
#define TMPN(x) ((tmpn_struct *) (x))

static int
tmpn_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    tmpn_ctx_struct * T = TMPN_CTX(ctx);
    int status = GR_SUCCESS;
    status |= gr_stream_write(out, T->is_signed ?
        "Transformed signed integers (fft_small, bits_bound " :
        "Transformed unsigned integers (fft_small, bits_bound ");
    status |= gr_stream_write_si(out, T->bits_bound);
    status |= gr_stream_write(out, ", terms_bound ");
    status |= gr_stream_write_ui(out, T->terms_bound);
    status |= gr_stream_write(out, ")");
    return status;
}

static void
tmpn_ctx_clear(gr_ctx_t ctx)
{
    tmpn_ctx_struct * T = TMPN_CTX(ctx);
    fft_small_plan_clear(T->P);
    fft_small_op_pool_clear(T->pool);
    flint_free(T);
}

static void
tmpn_init(gr_ptr x, gr_ctx_t ctx)
{
    tmpn_ctx_struct * T = TMPN_CTX(ctx);
    double * d = fft_small_op_pool_take(T->pool);

    if (d != NULL)
        fft_small_op_init_borrowed(&TMPN(x)->op, T->P, d);
    else
        fft_small_op_init(&TMPN(x)->op, T->P);
    TMPN(x)->nchunks = 0;
    TMPN(x)->small = 0;
    TMPN(x)->terms = 0;
    TMPN(x)->depth = 0;
    TMPN(x)->sign = 0;
    TMPN(x)->negs = 0;
}

/* element on caller-provided storage: 'data' must hold
   fft_small_op_sizeof_data of the context's plan, 4096-aligned, and
   outlive the element; gr_clear will not free it */

/* bytes of storage one element of this context needs */
ulong
gr_transformed_mpn_sizeof_data(gr_ctx_t ctx)
{
    return fft_small_op_sizeof_data(TMPN_CTX(ctx)->P);
}

static void
tmpn_clear(gr_ptr x, gr_ctx_t ctx)
{
    tmpn_ctx_struct * T = TMPN_CTX(ctx);
    char * d = (char *) TMPN(x)->op.data;

    /* a repeat clear is inert. Callers that release elements early,
       to free a slab for conversion scratch, clear them again through
       the vector clear at the end; returning a slab twice would hand
       the same storage to two later elements. */
    if (d == NULL)
        return;

    if (fft_small_op_pool_owns(T->pool, d))
    {
        fft_small_op_pool_give(T->pool, (double *) d);
        TMPN(x)->op.data = NULL;
        return;
    }
    fft_small_op_clear(&TMPN(x)->op);
    TMPN(x)->op.data = NULL;
    TMPN(x)->op.owns_data = 0;
}

static void
tmpn_swap(gr_ptr x, gr_ptr y, gr_ctx_t FLINT_UNUSED(ctx))
{
    FLINT_SWAP(tmpn_struct, *TMPN(x), *TMPN(y));
}

/* writes the represented value (reconstructed nondestructively);
   costs an inverse transform, which only diagnostics and the test
   framework pay */
static int
tmpn_write(gr_stream_t out, gr_srcptr x, gr_ctx_t ctx)
{
    slong need = gr_transformed_mpn_get_limbs(ctx, x);
    nn_ptr t;
    slong zn;
    int sgn, status;
    fmpz_t v;

    if (need <= 0)
    {
        return gr_stream_write(out, "0");
    }
    t = flint_malloc(need * sizeof(ulong));
    status = gr_transformed_mpn_get(t, need, &zn, &sgn, x, ctx);
    if (status == GR_SUCCESS)
    {
        fmpz_init(v);
        if (zn == 0)
            fmpz_zero(v);
        else
        {
            fmpz_set_ui_array(v, t, zn);
            if (sgn)
                fmpz_neg(v, v);
        }
        status |= gr_stream_write_fmpz(out, v);
        fmpz_clear(v);
    }
    flint_free(t);
    return status;
}

static int
tmpn_one(gr_ptr res, gr_ctx_t ctx)
{
    ulong v = 1;
    return gr_transformed_mpn_set(res, &v, 1, 0, ctx);
}

/* random depth-1 elements within the context's operand capacity;
   gr_test_ring drives the ring through these */
static int
tmpn_randtest(gr_ptr res, flint_rand_t state, gr_ctx_t ctx)
{
    tmpn_ctx_struct * T = TMPN_CTX(ctx);
    slong maxn = FLINT_MAX(1, (slong) (T->bits_bound / (2 * FLINT_BITS)));
    slong n = n_randint(state, maxn + 1);
    int sgn = T->is_signed ? (int) n_randint(state, 2) : 0;
    nn_ptr t;
    slong i;
    int status;

    if (n == 0 || n_randint(state, 8) == 0)
        return gr_zero(res, ctx);

    t = flint_malloc(n * sizeof(ulong));
    for (i = 0; i < n; i++)
        t[i] = n_randtest(state);
    t[n - 1] += (t[n - 1] == 0);
    status = gr_transformed_mpn_set(res, t, n, sgn, ctx);
    flint_free(t);
    return status;
}

static int
tmpn_zero(gr_ptr x, gr_ctx_t FLINT_UNUSED(ctx))
{
    TMPN(x)->nchunks = 0;
    TMPN(x)->terms = 0;
    TMPN(x)->depth = 0;
    TMPN(x)->sign = 0;
    TMPN(x)->negs = 0;
    TMPN(x)->small = 0;
    return GR_SUCCESS;
}

static int
tmpn_set(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    tmpn_ctx_struct * T = TMPN_CTX(ctx);

    if (res == x)
        return GR_SUCCESS;

    if (TMPN(x)->nchunks > 0)
        memcpy(TMPN(res)->op.data, TMPN(x)->op.data,
               T->P->np * T->P->stride * sizeof(double));
    TMPN(res)->op.domain = TMPN(x)->op.domain;
    TMPN(res)->op.itrunc = TMPN(x)->op.itrunc;
    TMPN(res)->nchunks = TMPN(x)->nchunks;
    TMPN(res)->terms = TMPN(x)->terms;
    TMPN(res)->depth = TMPN(x)->depth;
    TMPN(res)->sign = TMPN(x)->sign;
    TMPN(res)->negs = TMPN(x)->negs;
    TMPN(res)->small = TMPN(x)->small;
    return GR_SUCCESS;
}

static truth_t
tmpn_is_zero(gr_srcptr x, gr_ctx_t FLINT_UNUSED(ctx))
{
    if (TMPN(x)->nchunks == 0)
        return (TMPN(x)->small == 0) ? T_TRUE : T_FALSE;
    return T_UNKNOWN;
}

static truth_t
tmpn_equal(gr_srcptr x, gr_srcptr y, gr_ctx_t FLINT_UNUSED(ctx))
{
    if (TMPN(x)->nchunks == 0 && TMPN(y)->nchunks == 0)
        return (TMPN(x)->small == TMPN(y)->small) ? T_TRUE : T_FALSE;
    return T_UNKNOWN;
}

static int
tmpn_neg(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    int status = tmpn_set(res, x, ctx);
    if (status != GR_SUCCESS)
        return status;
    if (!TMPN_CTX(ctx)->is_signed
        && (TMPN(res)->nchunks > 0 || TMPN(res)->small > 0))
        return GR_UNABLE;
    if (TMPN(res)->nchunks > 0)
        TMPN(res)->sign ^= 1;
    TMPN(res)->small = -TMPN(res)->small;
    return status;
}

/* conversion in: the packing never sees signs */
int
gr_transformed_mpn_set(gr_ptr res, nn_srcptr a, slong an, int sign,
                       gr_ctx_t ctx)
{
    tmpn_ctx_struct * T = TMPN_CTX(ctx);

    while (an > 0 && a[an - 1] == 0)
        an--;

    if (an == 0)
        return tmpn_zero(res, ctx);

    if (sign && !T->is_signed)
        return GR_DOMAIN;

    if (an == 1 && a[0] == 1)
    {
        /* the coefficient stays in the small side: conversion costs
           nothing, and pointwise operations against it degenerate to
           additions or integer bookkeeping */
        tmpn_zero(res, ctx);
        TMPN(res)->small = sign ? -1 : 1;
        TMPN(res)->terms = 1;
        return GR_SUCCESS;
    }
    if ((ulong) an * FLINT_BITS > (ulong) T->zcap * T->P->bits)
        return GR_DOMAIN;

    fft_small_fft_mpn(&TMPN(res)->op, a, (ulong) an, T->P);
    TMPN(res)->nchunks = (slong) n_cdiv((ulong) an * FLINT_BITS, T->P->bits);
    TMPN(res)->terms = 1;
    TMPN(res)->depth = 1;
    TMPN(res)->sign = sign ? 1 : 0;
    TMPN(res)->negs = 0;
    TMPN(res)->small = 0;
    return GR_SUCCESS;
}

static int
tmpn_set_fmpz(gr_ptr elem, const fmpz_t a, gr_ctx_t tctx)
{
    fmpz c = *a;

    if (!COEFF_IS_MPZ(c))
    {
        ulong v = FLINT_ABS(c);
        return gr_transformed_mpn_set(elem, &v, c != 0, c < 0, tctx);
    }
    else
    {
        mpz_srcptr m = COEFF_TO_PTR(c);
        return gr_transformed_mpn_set(elem, m->_mp_d, FLINT_ABS(m->_mp_size), m->_mp_size < 0, tctx);
    }
}

/* limbs the reconstruction of an element with m chunks writes */
static slong
_tmpn_export_limbs(const tmpn_ctx_struct * T, slong m, int negs)
{
    ulong nslots = n_min((ulong) m, T->P->ztrunc);

    if (negs)
        return (slong) fft_small_export_mpn_signed_limbs(T->P, nslots);
    return (slong) fft_small_export_mpn_limbs(T->P, T->terms_bound, nslots);
}

/* the chunk slots of an element alias in a cyclic plan when it has more
   chunks than the transform length; only windows inside the context's
   window are then exact */
static int
_tmpn_aliased(const tmpn_ctx_struct * T, const tmpn_struct * x)
{
    return x->nchunks > (slong) T->P->ztrunc;
}

/* the element's transform data in a writable op: its own when
   destroying, else a copy (freed by _tmpn_op_done) */
static void
_tmpn_op_work(fft_small_op_struct * tmp, const tmpn_struct * x,
              int destroy, const fft_small_plan_struct * P)
{
    *tmp = x->op;
    if (!destroy)
    {
        tmp->data = flint_aligned_alloc(FLINT_FFT_SMALL_ALIGNMENT,
                    fft_small_op_sizeof_data(P));
        memcpy(tmp->data, x->op.data, P->np * P->stride * sizeof(double));
    }
}

static void
_tmpn_op_done(fft_small_op_struct * tmp, const tmpn_struct * x, int destroy)
{
    if (!destroy && tmp->data != x->op.data)
        flint_aligned_free(tmp->data);
}

/* fold the small side value into the normalized magnitude (zout, *len)
   with sign *sign; zout has room for the carry */
static void
_tmpn_fold_small(nn_ptr zout, slong * len, int * sign, slong sv)
{
    ulong av = (ulong) FLINT_ABS(sv);
    int ssign = sv < 0;
    slong i = *len;

    if (sv == 0)
        return;

    if (i == 0)
    {
        zout[0] = av;
        i = 1;
        *sign = ssign;
    }
    else if (*sign == ssign)
    {
        /* the export bound's terms headroom guarantees room */
        ulong cy = mpn_add_1(zout, zout, i, av);
        if (cy != 0)
            zout[i++] = cy;
    }
    else if (i > 1 || zout[0] > av)
    {
        mpn_sub_1(zout, zout, i, av);
        while (i > 0 && zout[i - 1] == 0)
            i--;
    }
    else
    {
        /* |T| <= |small|: the sign crosses */
        zout[0] = av - zout[0];
        i = (zout[0] != 0);
        *sign = ssign;
    }

    *len = i;
}

/*
    Conversions out. All share the fft_small machinery of the complex
    products: the truncated inverse transform with the normalizers this
    ring defers (fft_small_ifft_trunc), the exports, and the staging of
    reconstructions wider than the caller's window (fft_small_stage_*,
    through the context's pool). With 'destroy' the element's own
    transform buffer is consumed in place rather than copied; copying
    one is a full pass over np * 2^depth doubles, a substantial share of
    the conversion cost, so callers finished with the element should use
    the destructive forms (the element must then only be cleared).
*/

/* z receives the magnitude (zn limbs), *zn_out its normalized length
   and *sign the sign; a window short of the reconstruction is staged */
static int
_tmpn_get(nn_ptr z, slong zn, slong * zn_out, int * sign,
          gr_srcptr x, int destroy, gr_ctx_t ctx)
{
    tmpn_ctx_struct * T = TMPN_CTX(ctx);
    const fft_small_plan_struct * P = T->P;
    slong m = TMPN(x)->nchunks;
    int negs = TMPN(x)->negs;
    fft_small_op_struct tmp;
    nn_ptr zout;
    slong need, zoutn, i;

    if (m == 0)
    {
        slong sv = TMPN(x)->small;
        if (sv == 0)
        {
            if (zn > 0)
                flint_mpn_zero(z, zn);
            *zn_out = 0;
            *sign = 0;
            return GR_SUCCESS;
        }
        if (zn < 1)
            return GR_DOMAIN;
        z[0] = (ulong) FLINT_ABS(sv);
        if (zn > 1)
            flint_mpn_zero(z + 1, zn - 1);
        *zn_out = 1;
        *sign = sv < 0;
        return GR_SUCCESS;
    }

    if (_tmpn_aliased(T, TMPN(x)))
        return GR_UNABLE;

    need = _tmpn_export_limbs(T, m, negs);
    zout = fft_small_stage_begin(z, zn, need, T->pool, T->scratch_from_slab);
    zoutn = FLINT_MAX(zn, need);

    _tmpn_op_work(&tmp, TMPN(x), destroy, P);
    /* the unsigned reconstruction reads whole chunks up to the limb
       count, so the slots past the transformed range are cleared (they
       stand for zero chunks); the signed export reads m slots */
    fft_small_ifft_trunc(&tmp, (ulong) m,
            negs ? 0 : n_cdiv((ulong) need * FLINT_BITS, P->bits),
            T->m_orig, P);

    if (negs)
    {
        /* per-slot centered lifts: signed chunk values convert directly */
        int esign;
        fft_small_export_mpn_signed(zout, (ulong) need, &esign, &tmp,
                                    (ulong) m, P);
        *sign = TMPN(x)->sign ^ esign;
    }
    else
    {
        fft_small_export_mpn(zout, (ulong) need, &tmp, P);
        *sign = TMPN(x)->sign;
    }
    _tmpn_op_done(&tmp, TMPN(x), destroy);

    if (need < zoutn)
        flint_mpn_zero(zout + need, zoutn - need);

    i = need;
    while (i > 0 && zout[i - 1] == 0)
        i--;

    /* the element represents (-1)^sign * T + small */
    _tmpn_fold_small(zout, &i, sign, TMPN(x)->small);

    if (!fft_small_stage_end(z, zn, zout, i, T->pool))
        return GR_DOMAIN;   /* undersized for the value itself */

    *zn_out = i;
    if (i == 0)
        *sign = 0;
    return GR_SUCCESS;
}

slong
gr_transformed_mpn_get_limbs(gr_ctx_t ctx, gr_srcptr x)
{
    return _tmpn_export_limbs(TMPN_CTX(ctx), TMPN(x)->nchunks,
                              TMPN(x)->negs);
}

/* an upper bound over every element of the ring: what
   gr_transformed_mpn_get_limbs can return at the context's full chunk
   capacity -- callers size conversion staging from this instead of
   reconstructing the representation's limb requirements themselves */
slong
gr_transformed_mpn_get_limbs_bound(gr_ctx_t ctx)
{
    return _tmpn_export_limbs(TMPN_CTX(ctx), TMPN_CTX(ctx)->zcap, 1);
}

/* limbs required for a truncated conversion returning the limbs of the
   value starting at limb lo */
slong
gr_transformed_mpn_get_limbs_trunc(gr_ctx_t ctx, gr_srcptr x, slong lo)
{
    slong needf = _tmpn_export_limbs(TMPN_CTX(ctx), TMPN(x)->nchunks, 1);
    return FLINT_MAX(needf - lo, 0);
}

/* Truncated conversion: z receives zn limbs of |value| starting at limb
   lo (zero padded at the top), *sign the sign, *zn_out the significant
   limb count within the window. The discarded low tail perturbs the
   lowest returned limb by at most one unit (the mulhigh contract); with
   lo = 0 the conversion is exact. */
static int
_tmpn_get_trunc(nn_ptr z, slong zn, slong * zn_out, int * sign,
                slong lo, gr_srcptr x, int destroy, gr_ctx_t ctx)
{
    tmpn_ctx_struct * T = TMPN_CTX(ctx);
    const fft_small_plan_struct * P = T->P;
    slong m = TMPN(x)->nchunks;
    slong needf = _tmpn_export_limbs(T, m, 1);
    fft_small_op_struct tmp;
    nn_ptr zout;
    slong zoutn, i;
    int esign;

    if (m == 0 || lo >= needf)
    {
        if (m == 0 && lo == 0)
            return _tmpn_get(z, zn, zn_out, sign, x, destroy, ctx);
        flint_mpn_zero(z, zn);
        *zn_out = 0;
        *sign = 0;
        return GR_SUCCESS;
    }

    if (_tmpn_aliased(T, TMPN(x)))
        return GR_UNABLE;

    /* stage when the window does not reach the top of the
       reconstruction */
    zout = fft_small_stage_begin(z, zn, needf - lo, T->pool,
                                 T->scratch_from_slab);
    zoutn = FLINT_MAX(zn, needf - lo);

    _tmpn_op_work(&tmp, TMPN(x), destroy, P);
    fft_small_ifft_trunc(&tmp, (ulong) m, 0, T->m_orig, P);
    fft_small_export_mpn_signed_trunc(zout, (ulong) zoutn, &esign, &tmp,
                                      (ulong) m, (ulong) lo, P);
    _tmpn_op_done(&tmp, TMPN(x), destroy);

    *sign = TMPN(x)->sign ^ esign;

    i = zoutn;
    while (i > 0 && zout[i - 1] == 0)
        i--;

    /* for lo = 0 the window is all of T and the fold must be exact; for
       lo > 0 it is dropped: |small| < 2^63 <= 2^(FLINT_BITS lo), so it
       moves the window by at most one unit -- inside the documented
       truncation error -- and cannot flip the sign of a nonzero window */
    if (lo == 0)
        _tmpn_fold_small(zout, &i, sign, TMPN(x)->small);

    if (!fft_small_stage_end(z, zn, zout, i, T->pool))
        return GR_DOMAIN;

    *zn_out = i;
    if (i == 0)
        *sign = 0;
    return GR_SUCCESS;
}

/* Residue conversion: z receives the value mod 2^(FLINT_BITS zn) in
   two's complement, exactly, whatever the size of the value. Only the
   chunks below limb zn are reconstructed: the unsigned export reduces
   mod 2^(FLINT_BITS zn) by construction, and the signed one reads only
   the slots starting below limb zn. */
static int
_tmpn_get_low(nn_ptr z, slong zn, gr_srcptr x, int destroy, gr_ctx_t ctx)
{
    tmpn_ctx_struct * T = TMPN_CTX(ctx);
    const fft_small_plan_struct * P = T->P;
    slong m = TMPN(x)->nchunks;
    slong sv = TMPN(x)->small;
    fft_small_op_struct tmp;

    if (zn <= 0)
        return GR_SUCCESS;

    if (m == 0)
    {
        flint_mpn_zero(z, zn);
    }
    else
    {
        if (_tmpn_aliased(T, TMPN(x)))
            return GR_UNABLE;

        _tmpn_op_work(&tmp, TMPN(x), destroy, P);
        fft_small_ifft_trunc(&tmp, (ulong) m,
                TMPN(x)->negs ? 0 : n_cdiv((ulong) zn * FLINT_BITS, P->bits),
                T->m_orig, P);

        if (TMPN(x)->negs)
            fft_small_export_mpn_signed_low(z, (ulong) zn, &tmp, (ulong) m, P);
        else
        {
            /* the unsigned export reads whole chunks up to zn limbs,
               at most the transform length */
            ulong en = FLINT_MIN((ulong) zn, (ulong) _tmpn_export_limbs(T, m, 0));
            fft_small_export_mpn(z, en, &tmp, P);
            if ((ulong) zn > en)
                flint_mpn_zero(z + en, zn - en);
        }
        _tmpn_op_done(&tmp, TMPN(x), destroy);

        if (TMPN(x)->sign)
            mpn_neg(z, z, zn);
    }

    /* the element is (-1)^sign T + small */
    if (sv > 0)
        mpn_add_1(z, z, zn, (ulong) sv);
    else if (sv < 0)
        mpn_sub_1(z, z, zn, -(ulong) sv);

    return GR_SUCCESS;
}

/* Window conversion: z receives the limbs [lo, lo + zn) of the value as
   a two's-complement window, within one unit of floor(value / B^lo)
   (two with a nonzero small side value, which is folded exactly only
   for lo = 0) and exact for lo = 0; the value need not fit. Signed
   contexts only. Elements whose slots alias in a cyclic plan convert
   only inside the context's window. */
static int
_tmpn_get_window(nn_ptr z, slong zn, slong lo, gr_srcptr x, int destroy,
                 gr_ctx_t ctx)
{
    tmpn_ctx_struct * T = TMPN_CTX(ctx);
    const fft_small_plan_struct * P = T->P;
    slong m = TMPN(x)->nchunks;
    slong sv = TMPN(x)->small;
    fft_small_op_struct tmp;

    if (zn <= 0)
        return GR_SUCCESS;

    if (!T->is_signed || lo < 0)
        return GR_UNABLE;

    if (m == 0)
        flint_mpn_zero(z, zn);
    else
    {
        ulong nslots = n_min((ulong) m, P->ztrunc);

        if (_tmpn_aliased(T, TMPN(x)) &&
            !(T->windowed && lo >= T->win_lo && lo + zn <= T->win_hi))
            return GR_UNABLE;

        _tmpn_op_work(&tmp, TMPN(x), destroy, P);
        /* the sign goes in before the reconstruction: negating a floor
           window afterwards would cost another unit */
        if (TMPN(x)->sign)
            fft_small_op_neg(&tmp, &tmp, P);
        fft_small_ifft_trunc(&tmp, nslots, 0, T->m_orig, P);
        fft_small_export_mpn_signed_window(z, (ulong) zn, &tmp, nslots,
                                           (ulong) lo, P);
        _tmpn_op_done(&tmp, TMPN(x), destroy);
    }

    if (lo == 0)
    {
        if (sv > 0)
            mpn_add_1(z, z, zn, (ulong) sv);
        else if (sv < 0)
            mpn_sub_1(z, z, zn, -(ulong) sv);
    }

    return GR_SUCCESS;
}

static int
tmpn_add(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx);

static int
_tmpn_addsub(gr_ptr res, gr_srcptr x, gr_srcptr y, int ysign_flip,
             gr_ctx_t ctx)
{
    tmpn_ctx_struct * T = TMPN_CTX(ctx);
    int sy = TMPN(y)->sign ^ ysign_flip;
    ulong terms;

    if (TMPN(y)->nchunks == 0)
    {
        /* y is a pure small value */
        slong sy_small = ysign_flip ? -TMPN(y)->small : TMPN(y)->small;
        ulong t = TMPN(x)->terms + TMPN(y)->terms;
        int status;
        if (t < TMPN(x)->terms || t > T->terms_bound)
            return GR_UNABLE;
        status = tmpn_set(res, x, ctx);
        if (status != GR_SUCCESS)
            return status;
        TMPN(res)->small += sy_small;
        TMPN(res)->terms = t;
        if (!T->is_signed && TMPN(res)->nchunks == 0
            && TMPN(res)->small < 0)
            return GR_UNABLE;
        return GR_SUCCESS;
    }
    if (TMPN(x)->nchunks == 0)
    {
        ulong t = TMPN(x)->terms + TMPN(y)->terms;
        int status;
        slong sx_small = TMPN(x)->small;
        if (t < TMPN(x)->terms || t > T->terms_bound)
            return GR_UNABLE;
        status = tmpn_set(res, y, ctx);
        if (status != GR_SUCCESS)
            return status;
        if (ysign_flip)
        {
            if (!T->is_signed)
                return GR_UNABLE;
            TMPN(res)->sign ^= 1;
            TMPN(res)->small = -TMPN(res)->small;
        }
        TMPN(res)->small += sx_small;
        TMPN(res)->terms = t;
        return GR_SUCCESS;
    }

    terms = TMPN(x)->terms + TMPN(y)->terms;
    if (terms < TMPN(x)->terms || terms > T->terms_bound)
        return GR_UNABLE;

    TMPN(x)->op.domain = TMPN(y)->op.domain = FFT_SMALL_OP_PRIMAL;
    if (TMPN(x)->sign == sy)
    {
        /* same signs: pointwise addition, sign preserved */
        fft_small_op_add(&TMPN(res)->op, &TMPN(x)->op, &TMPN(y)->op, T->P);
        TMPN(res)->negs = TMPN(x)->negs | TMPN(y)->negs;
    }
    else
    {
        /* opposite signs: pointwise subtraction; the chunk values may go
           negative and the sign is resolved at conversion out */
        if (!T->is_signed)
            return GR_UNABLE;
        fft_small_op_sub(&TMPN(res)->op, &TMPN(x)->op, &TMPN(y)->op, T->P);
        TMPN(res)->negs = 1;
    }
    TMPN(res)->op.domain = FFT_SMALL_OP_PRIMAL;
    TMPN(res)->sign = TMPN(x)->sign;
    TMPN(res)->nchunks = FLINT_MAX(TMPN(x)->nchunks, TMPN(y)->nchunks);
    TMPN(res)->terms = terms;
    TMPN(res)->depth = FLINT_MAX(TMPN(x)->depth, TMPN(y)->depth);
    TMPN(res)->small = TMPN(x)->small
        + (ysign_flip ? -TMPN(y)->small : TMPN(y)->small);
    return GR_SUCCESS;
}

static int
tmpn_add(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
{
    return _tmpn_addsub(res, x, y, 0, ctx);
}

static int
tmpn_sub(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
{
    return _tmpn_addsub(res, x, y, 1, ctx);
}

/* the per-product chunk-count factor is provisioned inside the plan (its
   bound criterion multiplies by the operand chunk count); element terms
   count accumulated products only */
static int
_tmpn_mul_terms(ulong * terms, const tmpn_struct * x, const tmpn_struct * y,
                ulong terms_bound)
{
    ulong hi, t;

    umul_ppmm(hi, t, x->terms, y->terms);
    if (hi != 0 || t > terms_bound)
        return 0;
    *terms = t;
    return 1;
}

static int
tmpn_mul(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
{
    tmpn_ctx_struct * T = TMPN_CTX(ctx);
    slong m;
    ulong terms;

    if (TMPN(x)->nchunks == 0 || TMPN(y)->nchunks == 0)
    {
        /* a pure small factor: zero annihilates; +-1 copies the other
           operand with a sign adjustment; larger smalls (sums of
           trivial coefficients) are left to the caller */
        gr_srcptr sm = (TMPN(x)->nchunks == 0) ? x : y;
        gr_srcptr ot = (sm == x) ? y : x;
        slong sv = TMPN(sm)->small;
        int status;

        if (sv == 0)
            return tmpn_zero(res, ctx);
        if (sv != 1 && sv != -1)
            return GR_UNABLE;
        if (!_tmpn_mul_terms(&terms, TMPN(x), TMPN(y), T->terms_bound))
            return GR_UNABLE;
        status = tmpn_set(res, ot, ctx);
        if (status != GR_SUCCESS)
            return status;
        if (sv == -1)
        {
            if (!T->is_signed)
                return GR_UNABLE;
            if (TMPN(res)->nchunks > 0)
                TMPN(res)->sign ^= 1;
            TMPN(res)->small = -TMPN(res)->small;
        }
        TMPN(res)->terms = terms;
        return GR_SUCCESS;
    }

    /* a transformed operand carrying a small side value cannot enter a
       pointwise product without materializing it; the drivers only
       multiply freshly converted operands, whose small side is zero */
    if (TMPN(x)->small != 0 || TMPN(y)->small != 0)
        return GR_UNABLE;

    m = TMPN(x)->nchunks + TMPN(y)->nchunks - 1;
    if (m > T->zcap ||
        !_tmpn_mul_terms(&terms, TMPN(x), TMPN(y), T->terms_bound))
        return GR_UNABLE;
    if (TMPN(x)->depth + TMPN(y)->depth > T->max_depth)
        return GR_UNABLE;

    TMPN(x)->op.domain = TMPN(y)->op.domain = FFT_SMALL_OP_PRIMAL;
    /* a square is the same pointwise pass over one operand stream
       instead of two; the bookkeeping below is unchanged, since every
       field it derives from x and y coincides */
    if (x == y)
        fft_small_op_sqr(&TMPN(res)->op, &TMPN(x)->op, T->P);
    else
        fft_small_op_mul(&TMPN(res)->op, &TMPN(x)->op, &TMPN(y)->op, T->P);
    TMPN(res)->op.domain = FFT_SMALL_OP_PRIMAL;
    TMPN(res)->nchunks = m;
    TMPN(res)->terms = terms;
    TMPN(res)->depth = TMPN(x)->depth + TMPN(y)->depth;
    TMPN(res)->sign = TMPN(x)->sign ^ TMPN(y)->sign;
    TMPN(res)->negs = TMPN(x)->negs | TMPN(y)->negs;
    TMPN(res)->small = 0;
    return GR_SUCCESS;
}

static int
_tmpn_addsubmul(gr_ptr res, gr_srcptr x, gr_srcptr y, int subflip,
                gr_ctx_t ctx)
{
    tmpn_ctx_struct * T = TMPN_CTX(ctx);
    int psign;
    slong m;
    ulong terms, t2;

    if ((TMPN(x)->nchunks == 0 && TMPN(x)->small == 0)
        || (TMPN(y)->nchunks == 0 && TMPN(y)->small == 0))
        return GR_SUCCESS;   /* a zero factor: nothing to accumulate */

    if (TMPN(x)->nchunks == 0 || TMPN(y)->nchunks == 0)
    {
        /* a pure small +-1 factor: the accumulated product is the
           other operand up to sign, so the pointwise multiplication
           becomes an addition or subtraction (or pure integer
           bookkeeping when the other operand is small too) */
        gr_srcptr sm = (TMPN(x)->nchunks == 0) ? x : y;
        gr_srcptr ot = (sm == x) ? y : x;
        slong sv = TMPN(sm)->small;
        int psub;

        if (sv != 1 && sv != -1)
            return GR_UNABLE;
        psub = subflip ^ (sv < 0);
        if (TMPN(ot)->nchunks == 0)
        {
            ulong t = TMPN(res)->terms + TMPN(ot)->terms;
            if (t < TMPN(res)->terms || t > T->terms_bound)
                return GR_UNABLE;
            TMPN(res)->small += psub ? -TMPN(ot)->small
                                     : TMPN(ot)->small;
            TMPN(res)->terms = t;
            if (!T->is_signed && TMPN(res)->nchunks == 0
                && TMPN(res)->small < 0)
                return GR_UNABLE;
            return GR_SUCCESS;
        }
        return psub ? tmpn_sub(res, res, ot, ctx)
                    : tmpn_add(res, res, ot, ctx);
    }

    if (TMPN(x)->small != 0 || TMPN(y)->small != 0)
        return GR_UNABLE;

    /* with res aliasing an operand, the per-call domain bookkeeping
       below cannot mark the same struct as both input and accumulator;
       route through the ring's own guarded multiply and add */
    if (res == x || res == y)
    {
        gr_ptr t;
        int status;
        GR_TMP_INIT(t, ctx);
        status = tmpn_mul(t, x, y, ctx);
        if (status == GR_SUCCESS)
            status = subflip ? tmpn_sub(res, res, t, ctx)
                             : tmpn_add(res, res, t, ctx);
        GR_TMP_CLEAR(t, ctx);
        return status;
    }

    if (TMPN(res)->nchunks == 0)
    {
        /* the accumulator has no transform data, but it may hold an
           accumulated small side value, which must survive the
           product being written over it */
        slong rs = TMPN(res)->small;
        ulong rt = TMPN(res)->terms;
        int status = tmpn_mul(res, x, y, ctx);
        if (status != GR_SUCCESS)
            return status;
        if (subflip)
        {
            if (!T->is_signed)
                return GR_UNABLE;
            if (TMPN(res)->nchunks > 0)
                TMPN(res)->sign ^= 1;
            TMPN(res)->small = -TMPN(res)->small;
        }
        if (rs != 0)
        {
            ulong t = rt + TMPN(res)->terms;
            if (t < rt || t > T->terms_bound)
                return GR_UNABLE;
            TMPN(res)->small += rs;
            TMPN(res)->terms = t;
            if (!T->is_signed && TMPN(res)->nchunks == 0
                && TMPN(res)->small < 0)
                return GR_UNABLE;
        }
        return status;
    }

    m = TMPN(x)->nchunks + TMPN(y)->nchunks - 1;
    if (m > T->zcap ||
        !_tmpn_mul_terms(&t2, TMPN(x), TMPN(y), T->terms_bound))
        return GR_UNABLE;
    terms = TMPN(res)->terms + t2;
    if (terms < t2 || terms > T->terms_bound)
        return GR_UNABLE;
    if (TMPN(x)->depth + TMPN(y)->depth > T->max_depth)
        return GR_UNABLE;

    psign = TMPN(x)->sign ^ TMPN(y)->sign ^ subflip;

    TMPN(x)->op.domain = TMPN(y)->op.domain = FFT_SMALL_OP_PRIMAL;
    TMPN(res)->op.domain = FFT_SMALL_OP_PRODUCT;
    if (psign == TMPN(res)->sign)
    {
        /* accumulated product has the accumulator's sign: pointwise
           addmul */
        if (x == y)
            fft_small_op_addsqr(&TMPN(res)->op, &TMPN(x)->op, T->P);
        else
            fft_small_op_addmul(&TMPN(res)->op, &TMPN(x)->op, &TMPN(y)->op,
                                T->P);
        TMPN(res)->negs |= TMPN(x)->negs | TMPN(y)->negs;
    }
    else
    {
        /* opposite sign: the roles switch to a pointwise submul and the
           chunk values may go negative */
        if (!T->is_signed)
            return GR_UNABLE;
        if (x == y)
            fft_small_op_subsqr(&TMPN(res)->op, &TMPN(x)->op, T->P);
        else
            fft_small_op_submul(&TMPN(res)->op, &TMPN(x)->op, &TMPN(y)->op,
                                T->P);
        TMPN(res)->negs = 1;
    }
    TMPN(res)->op.domain = FFT_SMALL_OP_PRIMAL;
    TMPN(res)->nchunks = FLINT_MAX(TMPN(res)->nchunks, m);
    TMPN(res)->terms = terms;
    TMPN(res)->depth = FLINT_MAX(TMPN(res)->depth,
                                 TMPN(x)->depth + TMPN(y)->depth);
    return GR_SUCCESS;
}

static int
tmpn_addmul(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
{
    return _tmpn_addsubmul(res, x, y, 0, ctx);
}

static int
tmpn_submul(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
{
    return _tmpn_addsubmul(res, x, y, 1, ctx);
}

/* the generic square routes through gr_mul(res, x, x), which the
   coinciding-operand dispatch inside tmpn_mul already handles; this
   only removes the indirection and states the intent */
static int
tmpn_sqr(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    return tmpn_mul(res, x, x, ctx);
}

static gr_funcptr __tmpn_methods[GR_METHOD_TAB_SIZE];
static int __tmpn_methods_initialized = 0;

static gr_method_tab_input __tmpn_methods_input[] =
{
    {GR_METHOD_CTX_WRITE,       (gr_funcptr) (void (*)(void)) tmpn_ctx_write},
    {GR_METHOD_CTX_CLEAR,       (gr_funcptr) (void (*)(void)) tmpn_ctx_clear},
    {GR_METHOD_INIT,            (gr_funcptr) (void (*)(void)) tmpn_init},
    {GR_METHOD_CLEAR,           (gr_funcptr) (void (*)(void)) tmpn_clear},
    {GR_METHOD_SWAP,            (gr_funcptr) (void (*)(void)) tmpn_swap},
    {GR_METHOD_SET,             (gr_funcptr) (void (*)(void)) tmpn_set},
    {GR_METHOD_ZERO,            (gr_funcptr) (void (*)(void)) tmpn_zero},
    {GR_METHOD_WRITE,           (gr_funcptr) (void (*)(void)) tmpn_write},
    {GR_METHOD_ONE,             (gr_funcptr) (void (*)(void)) tmpn_one},
    {GR_METHOD_SET_FMPZ,        (gr_funcptr) (void (*)(void)) tmpn_set_fmpz},
    {GR_METHOD_RANDTEST,        (gr_funcptr) (void (*)(void)) tmpn_randtest},
    {GR_METHOD_IS_ZERO,         (gr_funcptr) (void (*)(void)) tmpn_is_zero},
    {GR_METHOD_EQUAL,           (gr_funcptr) (void (*)(void)) tmpn_equal},
    {GR_METHOD_NEG,             (gr_funcptr) (void (*)(void)) tmpn_neg},
    {GR_METHOD_ADD,             (gr_funcptr) (void (*)(void)) tmpn_add},
    {GR_METHOD_SUB,             (gr_funcptr) (void (*)(void)) tmpn_sub},
    {GR_METHOD_MUL,             (gr_funcptr) (void (*)(void)) tmpn_mul},
    {GR_METHOD_SQR,             (gr_funcptr) (void (*)(void)) tmpn_sqr},
    {GR_METHOD_ADDMUL,          (gr_funcptr) (void (*)(void)) tmpn_addmul},
    {GR_METHOD_SUBMUL,          (gr_funcptr) (void (*)(void)) tmpn_submul},
    {0,                         (gr_funcptr) (void (*)(void)) NULL},
};


int
gr_transformed_mpn_use_fit_buffer(gr_ctx_t ctx, slong num_live)
{
    tmpn_ctx_struct * T = TMPN_CTX(ctx);
    ulong stage;

    /* one reservation per context */
    if (T->pool->base != NULL)
        return GR_SUCCESS;
    if (num_live < 1)
        return GR_UNABLE;

    /* Conversion staging rides along in the same reservation, past the
       slabs. It is a fraction of one of them -- the reconstruction
       needs get_limbs_bound limbs against a slab's np * stride doubles
       -- so covering the ring's own conversions costs well under the
       extra element it would take from the pool. A caller that has
       promised a dead element instead (GR_TRANSFORMED_MPN_SCRATCH_FROM_SLAB)
       pays nothing here. Concurrent conversions beyond the one the
       region serves take plain allocations.

       The reserved head is zero: the ring's own operations request no
       scratch from the mpn_ctx, and interleaved requests from unrelated
       code on this thread are served from its secondary buffer. */
    stage = T->scratch_from_slab ? 0 :
        (ulong) gr_transformed_mpn_get_limbs_bound(ctx) * sizeof(ulong);

    if (!fft_small_op_pool_fit(T->pool, get_default_mpn_ctx(),
            fft_small_op_sizeof_data(T->P), num_live, stage,
            FFT_SMALL_OP_POOL_RESERVE))
        return GR_UNABLE;

    T->alloc_strategy = GR_TRANSFORMED_MPN_ALLOC_FIT_BUFFER;
    return GR_SUCCESS;
}

static int
_tmpn_ctx_init(gr_ctx_t ctx, slong bits_bound, slong terms_bound,
               int is_signed, slong num_live, int alloc_strategy,
               int windowed, slong lo, slong hi)
{
    tmpn_ctx_struct * T;
    ulong opn;

    if (num_live < 1)
        num_live = 4;

    if (bits_bound < FLINT_BITS || terms_bound < 1 ||
            (ulong) terms_bound > UWORD_MAX / 8)
        return GR_UNABLE;

    if (windowed && (lo < 0 || hi <= lo))
        return GR_UNABLE;

    T = FLINT_ARRAY_ALLOC(1, tmpn_ctx_struct);
    T->bits_bound = bits_bound;
    T->terms_bound = (ulong) terms_bound;
    T->max_depth = 2;
    T->is_signed = is_signed != 0;

    /* operands of up to bits_bound/2 bits each on either side (unbalanced
       operands simply occupy fewer chunks); provision the accumulation
       bound times two in the signed case for the bias headroom (biased
       chunk values reach twice the magnitude bound). The extra limbs
       give the plan chunk headroom for the top digits of the bias
       integer, whose per-chunk constant spans several chunks. A window
       plan (shared with the complex products) restricts the
       reconstruction to the chunks reaching [lo, hi) and may make the
       transforms cyclic. */
    opn = n_cdiv((ulong) bits_bound, 2 * FLINT_BITS) + 1;
    if (windowed ?
        !fft_small_plan_init_mpn_window(T->P, get_default_mpn_ctx(), opn, opn,
            (ulong) terms_bound, is_signed, (ulong) lo, (ulong) hi) :
        !fft_small_plan_init_mpn(T->P, get_default_mpn_ctx(), opn, opn,
            (ulong) terms_bound, is_signed))
    {
        flint_free(T);
        return GR_UNABLE;
    }

    /* decline when the declared number of simultaneously live elements
       would exceed the transform-storage budget, so callers fall back to
       slower algorithms instead of exhausting memory */
    if (fft_small_op_sizeof_data(T->P) >
            flint_fft_small_max_transformed_ring_size / (ulong) num_live)
    {
        fft_small_plan_clear(T->P);
        flint_free(T);
        return GR_UNABLE;
    }

    T->zcap = (slong) T->P->zn;
    T->windowed = windowed;
    T->win_lo = windowed ? lo : 0;
    T->win_hi = windowed ? hi : WORD_MAX;

    {
        ulong i;
        for (i = 0; i < T->P->np; i++)
        {
            T->m_orig[i] = T->P->m[i];
            T->P->m[i] = 1;
        }
    }

    /* Signed output conversion needs no per-context state: the signed
       export interprets each chunk slot's reconstruction as a centered
       residue, which is valid whenever chunk magnitudes stay below half
       the prime product -- guaranteed by the times-four provisioning of
       the accumulation bound above. */

    ctx->which_ring = GR_CTX_GR_TRANSFORMED_MPN;
    ctx->sizeof_elem = sizeof(tmpn_struct);
    ctx->size_limit = WORD_MAX;
    GR_CTX_DATA_AS_PTR(ctx) = T;
    ctx->methods = __tmpn_methods;

    /* concurrent initializations write identical data, which is the
       accepted pattern for gr method tables */
    if (!__tmpn_methods_initialized)
    {
        gr_method_tab_init(__tmpn_methods, __tmpn_methods_input);
        __tmpn_methods_initialized = 1;
    }

    T->alloc_strategy = alloc_strategy & GR_TRANSFORMED_MPN_ALLOC_STRATEGY_MASK;
    T->scratch_from_slab =
        (alloc_strategy & GR_TRANSFORMED_MPN_SCRATCH_FROM_SLAB) != 0;
    fft_small_op_pool_init(T->pool);
    if (T->alloc_strategy == GR_TRANSFORMED_MPN_ALLOC_FIT_BUFFER)
    {
        /* the operand slabs in the stable reserved tail; ring
           operations request no scratch from the mpn_ctx while the
           reservation is live */
        if (gr_transformed_mpn_use_fit_buffer(ctx, num_live)
                != GR_SUCCESS)
        {
            /* another reservation is live (a nested ring context):
               degrade to plain allocation */
            T->alloc_strategy = GR_TRANSFORMED_MPN_ALLOC_MALLOC;
        }
    }
    return GR_SUCCESS;
}

int
gr_ctx_init_transformed_mpn(gr_ctx_t ctx, slong bits_bound,
                            slong terms_bound, int is_signed, slong num_live,
                            int alloc_strategy)
{
    return _tmpn_ctx_init(ctx, bits_bound, terms_bound, is_signed, num_live,
                          alloc_strategy, 0, 0, 0);
}

int
gr_ctx_init_transformed_mpn_window(gr_ctx_t ctx, slong bits_bound,
                            slong terms_bound, int is_signed, slong num_live,
                            int alloc_strategy, slong lo, slong hi)
{
    return _tmpn_ctx_init(ctx, bits_bound, terms_bound, is_signed, num_live,
                          alloc_strategy, 1, lo, hi);
}


int
gr_transformed_mpn_get(nn_ptr z, slong zn, slong * zn_out, int * sign,
                       gr_srcptr x, gr_ctx_t ctx)
{
    return _tmpn_get(z, zn, zn_out, sign, x, 0, ctx);
}

/* as above but consumes x: only gr_clear may follow */
/* Convert x out into an fmpz, destructively. The width the
   reconstruction needs is the ring's own business, so it is taken from
   the element here and the destination grown to match: the conversion
   then writes straight into f's limbs, with no staging and no copy of
   the result. lo_limbs > 0 selects the truncated conversion. */
int
gr_transformed_mpn_get_fmpz_destructive(fmpz_t f, slong lo_limbs,
                                        gr_ptr x, gr_ctx_t ctx)
{
    slong w, zn;
    int sg, status;
    mpz_ptr m;

    w = (lo_limbs == 0) ? gr_transformed_mpn_get_limbs(ctx, x)
                        : gr_transformed_mpn_get_limbs_trunc(ctx, x, lo_limbs);
    if (w <= 0)
        w = 1;

    m = _fmpz_promote(f);
    if (lo_limbs == 0)
        status = gr_transformed_mpn_get_destructive(FLINT_MPZ_REALLOC(m, w),
                     w, &zn, &sg, x, ctx);
    else
        status = gr_transformed_mpn_get_trunc_destructive(
                     FLINT_MPZ_REALLOC(m, w), w, &zn, &sg, lo_limbs, x, ctx);

    if (status != GR_SUCCESS)
    {
        _fmpz_demote(f);
        fmpz_zero(f);
        return status;
    }

    m->_mp_size = sg ? -zn : zn;
    _fmpz_demote_val(f);
    return GR_SUCCESS;
}

int
gr_transformed_mpn_get_destructive(nn_ptr z, slong zn, slong * zn_out,
                                   int * sign, gr_ptr x, gr_ctx_t ctx)
{
    return _tmpn_get(z, zn, zn_out, sign, x, 1, ctx);
}

int
gr_transformed_mpn_get_trunc(nn_ptr z, slong zn, slong * zn_out, int * sign,
                             slong lo, gr_srcptr x, gr_ctx_t ctx)
{
    return _tmpn_get_trunc(z, zn, zn_out, sign, lo, x, 0, ctx);
}

int
gr_transformed_mpn_get_trunc_destructive(nn_ptr z, slong zn, slong * zn_out,
                                         int * sign, slong lo, gr_ptr x,
                                         gr_ctx_t ctx)
{
    return _tmpn_get_trunc(z, zn, zn_out, sign, lo, x, 1, ctx);
}

int
gr_transformed_mpn_get_low(nn_ptr z, slong zn, gr_srcptr x, gr_ctx_t ctx)
{
    return _tmpn_get_low(z, zn, x, 0, ctx);
}

int
gr_transformed_mpn_get_low_destructive(nn_ptr z, slong zn, gr_ptr x,
                                       gr_ctx_t ctx)
{
    return _tmpn_get_low(z, zn, x, 1, ctx);
}

int
gr_transformed_mpn_get_window(nn_ptr z, slong zn, slong lo, gr_srcptr x,
                              gr_ctx_t ctx)
{
    return _tmpn_get_window(z, zn, lo, x, 0, ctx);
}

int
gr_transformed_mpn_get_window_destructive(nn_ptr z, slong zn, slong lo,
                                          gr_ptr x, gr_ctx_t ctx)
{
    return _tmpn_get_window(z, zn, lo, x, 1, ctx);
}

#else /* FLINT_HAVE_FFT_SMALL */

int
gr_ctx_init_transformed_mpn_window(gr_ctx_t FLINT_UNUSED(ctx), slong FLINT_UNUSED(bits_bound),
                            slong FLINT_UNUSED(terms_bound), int FLINT_UNUSED(is_signed), slong FLINT_UNUSED(num_live),
                            int FLINT_UNUSED(alloc_strategy), slong FLINT_UNUSED(lo), slong FLINT_UNUSED(hi))
{
    return GR_UNABLE;
}

int
gr_transformed_mpn_get_window(nn_ptr FLINT_UNUSED(z), slong FLINT_UNUSED(zn), slong FLINT_UNUSED(lo),
                              gr_srcptr FLINT_UNUSED(x), gr_ctx_t FLINT_UNUSED(ctx))
{
    return GR_UNABLE;
}

int
gr_transformed_mpn_get_window_destructive(nn_ptr FLINT_UNUSED(z), slong FLINT_UNUSED(zn), slong FLINT_UNUSED(lo),
                                          gr_ptr FLINT_UNUSED(x), gr_ctx_t FLINT_UNUSED(ctx))
{
    return GR_UNABLE;
}

/* Without fft_small the transformed representation does not exist: the
   constructor reports GR_UNABLE and callers fall through to their plain
   code paths. The remaining functions are unreachable then, but the
   full public surface is defined so the library links; each returns the
   neutral answer for its type. */

int
gr_ctx_init_transformed_mpn(gr_ctx_t FLINT_UNUSED(ctx), slong FLINT_UNUSED(bits_bound),
                            slong FLINT_UNUSED(terms_bound), int FLINT_UNUSED(is_signed), slong FLINT_UNUSED(num_live), int FLINT_UNUSED(alloc_strategy))
{
    return GR_UNABLE;
}

int
gr_transformed_mpn_set(gr_ptr FLINT_UNUSED(res), nn_srcptr FLINT_UNUSED(a), slong FLINT_UNUSED(an), int FLINT_UNUSED(sign),
                       gr_ctx_t FLINT_UNUSED(ctx))
{
    return GR_UNABLE;
}

int
gr_transformed_mpn_get(nn_ptr FLINT_UNUSED(z), slong FLINT_UNUSED(zn), slong * FLINT_UNUSED(zn_out), int * FLINT_UNUSED(sign),
                       gr_srcptr FLINT_UNUSED(x), gr_ctx_t FLINT_UNUSED(ctx))
{
    return GR_UNABLE;
}

int
gr_transformed_mpn_get_destructive(nn_ptr FLINT_UNUSED(z), slong FLINT_UNUSED(zn), slong * FLINT_UNUSED(zn_out),
                                   int * FLINT_UNUSED(sign), gr_ptr FLINT_UNUSED(x), gr_ctx_t FLINT_UNUSED(ctx))
{
    return GR_UNABLE;
}

int
gr_transformed_mpn_get_trunc(nn_ptr FLINT_UNUSED(z), slong FLINT_UNUSED(zn), slong * FLINT_UNUSED(zn_out), int * FLINT_UNUSED(sign),
                             slong FLINT_UNUSED(lo), gr_srcptr FLINT_UNUSED(x), gr_ctx_t FLINT_UNUSED(ctx))
{
    return GR_UNABLE;
}

int
gr_transformed_mpn_get_trunc_destructive(nn_ptr FLINT_UNUSED(z), slong FLINT_UNUSED(zn), slong * FLINT_UNUSED(zn_out),
                                         int * FLINT_UNUSED(sign), slong FLINT_UNUSED(lo), gr_ptr FLINT_UNUSED(x),
                                         gr_ctx_t FLINT_UNUSED(ctx))
{
    return GR_UNABLE;
}

int
gr_transformed_mpn_get_low(nn_ptr FLINT_UNUSED(z), slong FLINT_UNUSED(zn), gr_srcptr FLINT_UNUSED(x), gr_ctx_t FLINT_UNUSED(ctx))
{
    return GR_UNABLE;
}

int
gr_transformed_mpn_get_low_destructive(nn_ptr FLINT_UNUSED(z), slong FLINT_UNUSED(zn), gr_ptr FLINT_UNUSED(x), gr_ctx_t FLINT_UNUSED(ctx))
{
    return GR_UNABLE;
}

slong
gr_transformed_mpn_get_limbs(gr_ctx_t FLINT_UNUSED(ctx), gr_srcptr FLINT_UNUSED(x))
{
    return 0;
}

slong
gr_transformed_mpn_get_limbs_bound(gr_ctx_t FLINT_UNUSED(ctx))
{
    return 0;
}

slong
gr_transformed_mpn_get_limbs_trunc(gr_ctx_t FLINT_UNUSED(ctx), gr_srcptr FLINT_UNUSED(x), slong FLINT_UNUSED(lo))
{
    return 0;
}

int
gr_transformed_mpn_use_fit_buffer(gr_ctx_t FLINT_UNUSED(ctx), slong FLINT_UNUSED(num_live))
{
    return GR_UNABLE;
}

ulong
gr_transformed_mpn_sizeof_data(gr_ctx_t FLINT_UNUSED(ctx))
{
    return 0;
}

#endif
