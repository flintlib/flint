/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
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

    A context is read-only after construction and may be shared by any
    number of threads; a single element belongs to one thread at a time.
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
} tmpn_ctx_struct;

typedef struct
{
    fft_small_op_struct op;
    slong nchunks;                  /* virtual chunk length; 0 = zero */
    ulong terms;
    ulong depth;
    int sign;                       /* 0 or 1: (-1)^sign */
    int negs;                       /* chunk values may be signed */
} tmpn_struct;

#define TMPN_CTX(ctx) ((tmpn_ctx_struct *) GR_CTX_DATA_AS_PTR(ctx))
#define TMPN(x) ((tmpn_struct *) (x))

/* d = m*d + b over the first nblk blocks (see the polynomial rings) */
static void
_tmpn_scale(const sd_fft_ctx_struct * Q, double * d, ulong m_, ulong nblk)
{
    vec8d m = vec8d_set_d(vec1d_reduce_0n_to_pmhn((slong) m_, Q->p));
    vec8d n = vec8d_set_d(Q->p);
    vec8d ninv = vec8d_set_d(Q->pinv);
    ulong I;

    for (I = 0; I < nblk; I++)
    {
        double * dx = d + sd_fft_ctx_blk_offset(I);
        ulong j = 0; do {
            vec8d x0, x1;
            x0 = vec8d_load(dx + j + 0);
            x1 = vec8d_load(dx + j + 8);
            x0 = vec8d_mulmod(x0, m, n, ninv);
            x1 = vec8d_mulmod(x1, m, n, ninv);
            vec8d_store(dx + j + 0, x0);
            vec8d_store(dx + j + 8, x1);
        } while (j += 16, j < BLK_SZ);
    }
}

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
    flint_free(T);
}

static void
tmpn_init(gr_ptr x, gr_ctx_t ctx)
{
    fft_small_op_init(&TMPN(x)->op, TMPN_CTX(ctx)->P);
    TMPN(x)->nchunks = 0;
    TMPN(x)->terms = 0;
    TMPN(x)->depth = 0;
    TMPN(x)->sign = 0;
    TMPN(x)->negs = 0;
}

/* element on caller-provided storage: 'data' must hold
   fft_small_op_sizeof_data of the context's plan, 4096-aligned, and
   outlive the element; gr_clear will not free it */
void
gr_transformed_mpn_init_borrowed(gr_ptr x, double * data, gr_ctx_t ctx)
{
    fft_small_op_init_borrowed(&TMPN(x)->op, TMPN_CTX(ctx)->P, data);
    TMPN(x)->nchunks = 0;
    TMPN(x)->terms = 0;
    TMPN(x)->depth = 0;
    TMPN(x)->sign = 0;
    TMPN(x)->negs = 0;
}

/* bytes of storage one element of this context needs */
ulong
gr_transformed_mpn_sizeof_data(gr_ctx_t ctx)
{
    return fft_small_op_sizeof_data(TMPN_CTX(ctx)->P);
}

static void
tmpn_clear(gr_ptr x, gr_ctx_t FLINT_UNUSED(ctx))
{
    fft_small_op_clear(&TMPN(x)->op);
}

static void
tmpn_swap(gr_ptr x, gr_ptr y, gr_ctx_t FLINT_UNUSED(ctx))
{
    FLINT_SWAP(tmpn_struct, *TMPN(x), *TMPN(y));
}

static int
tmpn_zero(gr_ptr x, gr_ctx_t FLINT_UNUSED(ctx))
{
    TMPN(x)->nchunks = 0;
    TMPN(x)->terms = 0;
    TMPN(x)->depth = 0;
    TMPN(x)->sign = 0;
    TMPN(x)->negs = 0;
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
    return GR_SUCCESS;
}

static truth_t
tmpn_is_zero(gr_srcptr x, gr_ctx_t FLINT_UNUSED(ctx))
{
    return (TMPN(x)->nchunks == 0) ? T_TRUE : T_UNKNOWN;
}

static truth_t
tmpn_equal(gr_srcptr x, gr_srcptr y, gr_ctx_t FLINT_UNUSED(ctx))
{
    if (TMPN(x)->nchunks == 0 && TMPN(y)->nchunks == 0)
        return T_TRUE;
    return T_UNKNOWN;
}

static int
tmpn_neg(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    int status = tmpn_set(res, x, ctx);
    if (status == GR_SUCCESS && TMPN(res)->nchunks > 0)
    {
        if (!TMPN_CTX(ctx)->is_signed)
            return GR_UNABLE;
        TMPN(res)->sign ^= 1;
    }
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
    if ((ulong) an * FLINT_BITS > (ulong) T->zcap * T->P->bits)
        return GR_DOMAIN;

    fft_small_fft_mpn(&TMPN(res)->op, a, (ulong) an, T->P);
    TMPN(res)->nchunks = (slong) n_cdiv((ulong) an * FLINT_BITS, T->P->bits);
    TMPN(res)->terms = 1;
    TMPN(res)->depth = 1;
    TMPN(res)->sign = sign ? 1 : 0;
    TMPN(res)->negs = 0;
    return GR_SUCCESS;
}

/* number of limbs certainly sufficient for the biased reconstruction of
   an element with m chunks */
static slong
_tmpn_export_limbs_raw(const fft_small_plan_struct * P, ulong terms_bound,
                       ulong mchunks)
{
    return (slong) n_cdiv(mchunks * P->bits + 2 * P->bits
            + FLINT_BIT_COUNT(terms_bound)
            + FLINT_BIT_COUNT((ulong) P->zn) + 4,
            FLINT_BITS);
}

static slong
_tmpn_export_limbs(const tmpn_ctx_struct * T, slong m, int negs)
{
    if (negs)
    {
        /* the signed export reads round_up(m) slots and requires
           64 zn >= nslots bits + 64 coeff_len + 2 */
        ulong clen = (T->P->crts + T->P->np - 1)->coeff_len;
        return (slong) n_cdiv((ulong) m * T->P->bits
                + FLINT_BITS * clen + 2 + (FLINT_BITS - 1), FLINT_BITS);
    }
    return _tmpn_export_limbs_raw(T->P, T->terms_bound,
            n_min(n_round_up((ulong) m, BLK_SZ), T->P->zn));
}

/* conversion out: z receives the magnitude (zn limbs allocated by the
   caller, at least gr_transformed_mpn_get_limbs(ctx, x)); *zn_out is the
   normalized limb count and *sign the sign bit */
/*
    Conversion out. With 'destroy' the element's own transform buffer is
    consumed in place rather than copied. Both the inverse transform and
    the reconstruction need a writable copy of the evaluations, and
    copying one is a full pass over np * 2^depth doubles -- 3 MB for a
    megabyte-size operand, a substantial share of the conversion cost.
    Callers finished with the element should use the destructive form;
    the element is then unusable and must only be cleared.
*/
static int
_tmpn_get(nn_ptr z, slong zn, slong * zn_out, int * sign,
          gr_srcptr x, int destroy, gr_ctx_t ctx)
{
    tmpn_ctx_struct * T = TMPN_CTX(ctx);
    const fft_small_plan_struct * P = T->P;
    slong m = TMPN(x)->nchunks;
    int negs = TMPN(x)->negs;
    slong need = _tmpn_export_limbs(T, m, negs);
    fft_small_op_struct tmp;
    double * sdata;
    ulong itr;
    slong i;

    if (m == 0)
    {
        *zn_out = 0;
        *sign = 0;
        return GR_SUCCESS;
    }

    if (zn < need)
        return GR_DOMAIN;

    if (destroy)
    {
        sdata = NULL;
        tmp = TMPN(x)->op;
    }
    else
    {
        sdata = flint_aligned_alloc(FLINT_FFT_SMALL_ALIGNMENT,
                n_round_up(P->np * P->stride * sizeof(double), FLINT_FFT_SMALL_ALIGNMENT));
        tmp = TMPN(x)->op;
        tmp.data = sdata;
        memcpy(sdata, TMPN(x)->op.data,
               P->np * P->stride * sizeof(double));
    }
    itr = n_round_up((ulong) m, BLK_SZ);

    for (i = 0; i < (slong) P->np; i++)
    {
        sd_fft_ctx_struct * Q = P->ffts + P->offset + i;
        double * d = tmp.data + P->stride * i;

        sd_ifft_trunc(Q, d, P->depth, itr);
        _tmpn_scale(Q, d, T->m_orig[i], itr / BLK_SZ);

        /* the unsigned reconstruction reads whole chunks up to the
           requested limb count; slots beyond the inverse transform's
           range are garbage and represent true zero chunks, so clear
           them (the signed export reads exactly itr slots instead) */
        if (!negs)
        {
            ulong chi = n_min(n_pow2(P->depth),
                              n_round_up(n_cdiv((ulong) need * FLINT_BITS,
                                                P->bits), BLK_SZ));
            ulong I;
            for (I = itr / BLK_SZ; I < chi / BLK_SZ; I++)
                memset(d + sd_fft_ctx_blk_offset(I), 0,
                       BLK_SZ * sizeof(double));
        }
    }
    tmp.domain = FFT_SMALL_OP_PRODUCT;
    if (negs)
    {
        /* per-slot centered lifts: signed chunk values convert directly,
           with no bias operand and no precomputation */
        int esign;
        fft_small_export_mpn_signed(z, (ulong) need, &esign, &tmp,
                                    (ulong) m, P);
        *sign = TMPN(x)->sign ^ esign;
    }
    else
    {
        fft_small_export_mpn(z, (ulong) need, &tmp, P);
        *sign = TMPN(x)->sign;
    }
    if (sdata != NULL)
        flint_aligned_free(sdata);

    if (need < zn)
        flint_mpn_zero(z + need, zn - need);

    i = need;
    while (i > 0 && z[i - 1] == 0)
        i--;
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
    double * sdata;
    ulong itr;
    slong i;
    int esign;

    if (m == 0 || lo >= needf)
    {
        flint_mpn_zero(z, zn);
        *zn_out = 0;
        *sign = 0;
        return GR_SUCCESS;
    }

    if (lo + zn < needf)
        return GR_DOMAIN;

    itr = n_round_up((ulong) m, BLK_SZ);

    if (destroy)
    {
        sdata = NULL;
        tmp = TMPN(x)->op;
    }
    else
    {
        sdata = flint_aligned_alloc(FLINT_FFT_SMALL_ALIGNMENT,
                n_round_up(P->np * P->stride * sizeof(double), FLINT_FFT_SMALL_ALIGNMENT));
        tmp = TMPN(x)->op;
        tmp.data = sdata;
        memcpy(sdata, TMPN(x)->op.data,
               P->np * P->stride * sizeof(double));
    }

    for (i = 0; i < (slong) P->np; i++)
    {
        sd_fft_ctx_struct * Q = P->ffts + P->offset + i;
        double * d = tmp.data + P->stride * i;

        sd_ifft_trunc(Q, d, P->depth, itr);
        _tmpn_scale(Q, d, T->m_orig[i], itr / BLK_SZ);
    }
    tmp.domain = FFT_SMALL_OP_PRODUCT;
    fft_small_export_mpn_signed_trunc(z, (ulong) zn, &esign, &tmp,
                                      (ulong) m, (ulong) lo, P);
    if (sdata != NULL)
        flint_aligned_free(sdata);

    *sign = TMPN(x)->sign ^ esign;

    i = zn;
    while (i > 0 && z[i - 1] == 0)
        i--;
    *zn_out = i;
    if (i == 0)
        *sign = 0;
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
        return tmpn_set(res, x, ctx);
    if (TMPN(x)->nchunks == 0)
    {
        int status = tmpn_set(res, y, ctx);
        if (status == GR_SUCCESS && (TMPN(res)->sign ^ ysign_flip) !=
                TMPN(res)->sign)
        {
            if (!T->is_signed)
                return GR_UNABLE;
            TMPN(res)->sign ^= ysign_flip;
        }
        return status;
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
        return tmpn_zero(res, ctx);

    m = TMPN(x)->nchunks + TMPN(y)->nchunks - 1;
    if (m > T->zcap ||
        !_tmpn_mul_terms(&terms, TMPN(x), TMPN(y), T->terms_bound))
        return GR_UNABLE;
    if (TMPN(x)->depth + TMPN(y)->depth > T->max_depth)
        return GR_UNABLE;

    TMPN(x)->op.domain = TMPN(y)->op.domain = FFT_SMALL_OP_PRIMAL;
    fft_small_op_mul(&TMPN(res)->op, &TMPN(x)->op, &TMPN(y)->op, T->P);
    TMPN(res)->op.domain = FFT_SMALL_OP_PRIMAL;
    TMPN(res)->nchunks = m;
    TMPN(res)->terms = terms;
    TMPN(res)->depth = TMPN(x)->depth + TMPN(y)->depth;
    TMPN(res)->sign = TMPN(x)->sign ^ TMPN(y)->sign;
    TMPN(res)->negs = TMPN(x)->negs | TMPN(y)->negs;
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

    if (TMPN(x)->nchunks == 0 || TMPN(y)->nchunks == 0)
        return GR_SUCCESS;

    if (TMPN(res)->nchunks == 0)
    {
        int status = tmpn_mul(res, x, y, ctx);
        if (status == GR_SUCCESS && subflip && TMPN(res)->nchunks > 0)
        {
            if (!T->is_signed)
                return GR_UNABLE;
            TMPN(res)->sign ^= 1;
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
        fft_small_op_addmul(&TMPN(res)->op, &TMPN(x)->op, &TMPN(y)->op, T->P);
        TMPN(res)->negs |= TMPN(x)->negs | TMPN(y)->negs;
    }
    else
    {
        /* opposite sign: the roles switch to a pointwise submul and the
           chunk values may go negative */
        if (!T->is_signed)
            return GR_UNABLE;
        fft_small_op_submul(&TMPN(res)->op, &TMPN(x)->op, &TMPN(y)->op, T->P);
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
    {GR_METHOD_IS_ZERO,         (gr_funcptr) (void (*)(void)) tmpn_is_zero},
    {GR_METHOD_EQUAL,           (gr_funcptr) (void (*)(void)) tmpn_equal},
    {GR_METHOD_NEG,             (gr_funcptr) (void (*)(void)) tmpn_neg},
    {GR_METHOD_ADD,             (gr_funcptr) (void (*)(void)) tmpn_add},
    {GR_METHOD_SUB,             (gr_funcptr) (void (*)(void)) tmpn_sub},
    {GR_METHOD_MUL,             (gr_funcptr) (void (*)(void)) tmpn_mul},
    {GR_METHOD_ADDMUL,          (gr_funcptr) (void (*)(void)) tmpn_addmul},
    {GR_METHOD_SUBMUL,          (gr_funcptr) (void (*)(void)) tmpn_submul},
    {0,                         (gr_funcptr) (void (*)(void)) NULL},
};

int
gr_ctx_init_transformed_mpn(gr_ctx_t ctx, slong bits_bound,
                            slong terms_bound, int is_signed, slong num_live)
{
    tmpn_ctx_struct * T;
    ulong opn;

    if (num_live < 1)
        num_live = 4;

    if (bits_bound < FLINT_BITS || terms_bound < 1 ||
            (ulong) terms_bound > UWORD_MAX / 8)
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
       integer, whose per-chunk constant spans several chunks */
    opn = n_cdiv((ulong) bits_bound, 2 * FLINT_BITS) + 1;
    if (!fft_small_plan_init_mpn(T->P, get_default_mpn_ctx(), opn, opn,
            (is_signed ? 2 : 1) * (ulong) terms_bound))
    {
        flint_free(T);
        return GR_UNABLE;
    }

    /* decline when the declared number of simultaneously live elements
       would exceed the transform-storage budget, so callers fall back to
       slower algorithms instead of exhausting memory */
    {
        ulong per = n_round_up(T->P->np * T->P->stride * sizeof(double),
                               FLINT_FFT_SMALL_ALIGNMENT);
        if (per > flint_fft_small_max_transformed_ring_size / (ulong) num_live)
        {
            fft_small_plan_clear(T->P);
            flint_free(T);
            return GR_UNABLE;
        }
    }

    T->zcap = (slong) T->P->zn;

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

    return GR_SUCCESS;
}


int
gr_transformed_mpn_get(nn_ptr z, slong zn, slong * zn_out, int * sign,
                       gr_srcptr x, gr_ctx_t ctx)
{
    return _tmpn_get(z, zn, zn_out, sign, x, 0, ctx);
}

/* as above but consumes x: only gr_clear may follow */
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

#else /* FLINT_HAVE_FFT_SMALL */

/* Without fft_small the transformed representation does not exist: the
   constructor reports GR_UNABLE and callers fall through to their plain
   code paths. The remaining functions are unreachable then, but the
   full public surface is defined so the library links; each returns the
   neutral answer for its type. */

int
gr_ctx_init_transformed_mpn(gr_ctx_t FLINT_UNUSED(ctx), slong FLINT_UNUSED(bits_bound),
                            slong FLINT_UNUSED(terms_bound), int FLINT_UNUSED(is_signed), slong FLINT_UNUSED(num_live))
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

slong
gr_transformed_mpn_get_limbs(gr_ctx_t FLINT_UNUSED(ctx), gr_srcptr FLINT_UNUSED(x))
{
    return 0;
}

slong
gr_transformed_mpn_get_limbs_trunc(gr_ctx_t FLINT_UNUSED(ctx), gr_srcptr FLINT_UNUSED(x), slong FLINT_UNUSED(lo))
{
    return 0;
}

void
gr_transformed_mpn_init_borrowed(gr_ptr FLINT_UNUSED(x), double * FLINT_UNUSED(data), gr_ctx_t FLINT_UNUSED(ctx))
{
}

ulong
gr_transformed_mpn_sizeof_data(gr_ctx_t FLINT_UNUSED(ctx))
{
    return 0;
}

#endif
