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
#include "gr_poly.h"
#include "nmod.h"
#include "mpn_mod.h"

#if FLINT_HAVE_FFT_SMALL

#include "machine_vectors.h"
#include "fft_small.h"

/* truncation granularity follows the plan's depth: whole blocks at
   block depths, the full (short) transform length below one block */

/* the anticipated granularity for a length before any plan exists */
#define _len_trunc(x) \
    ((n_clog2(x) < LG_BLK_SZ) ? n_pow2(n_max((ulong) 4, n_clog2(x))) \
                              : n_round_up((x), BLK_SZ))
#define _op_trunc(x, Pln) \
    (((Pln)->depth < LG_BLK_SZ) ? n_pow2((Pln)->depth) \
                                : n_round_up((x), BLK_SZ))

/*
    Ring of mpn_mod polynomials of bounded length in fft_small transformed
    representation, following gr/nmod_transformed_poly.c: elements are
    plain per-prime evaluations (the plan's pointwise normalizers are
    overridden to 1 and the saved originals applied once at conversion
    out), carrying a virtual length, an elementary product counter, a
    multiplicative depth and a negativity flag. Conversions in and out go
    through the per-coefficient reductions and chinese remaindering of the
    mpn_mod fused driver (fft_small_fft_mpn_mod /
    fft_small_export_mpn_mod_range).

    Differences from the nmod ring:
    - the chinese remaindering reconstructs full-width integers and asserts
      them below half the prime product, so the plan is provisioned with an
      accumulation bound of 4 * terms_bound (one bit for the signed-value
      bias, one for the unsigned-half constraint);
    - there is no exact single-prime path (the modulus never matches an
      fft prime), so the multiplicative depth is always capped at 2 and
      the bias always provisioned;
    - the base context is retained by pointer for the modular reductions
      at conversion out, so it must outlive this context.

    A context is read-only after construction and may be shared by any
    number of threads (a single element still belongs to one thread at a
    time).
*/

typedef struct
{
    fft_small_plan_t P;
    gr_ctx_struct * base;           /* must outlive this context */
    slong nlimbs;
    slong N;                        /* length capacity */
    ulong terms_bound;              /* accumulation capacity */
    ulong max_depth;                /* multiplicative depth capacity */
    ulong m_orig[MPN_CTX_NCRTS];    /* saved export normalizers */
    ulong bias_res[MPN_CTX_NCRTS];  /* C * inv(cop_i) mod p_i */
    nn_ptr bias_mod_n;              /* C mod n, nlimbs */
} mtpoly_ctx_struct;

typedef struct
{
    fft_small_op_struct op;
    slong len;                      /* virtual length; 0 = zero element */
    ulong terms;                    /* elementary product count bound */
    ulong depth;                    /* multiplicative depth in base elements */
    int negs;                       /* true integer coeffs may be negative */
} mtpoly_struct;

#define MTPOLY_CTX(ctx) ((mtpoly_ctx_struct *) GR_CTX_DATA_AS_PTR(ctx))
#define MTPOLY(x) ((mtpoly_struct *) (x))

/* d = m*d + b over the element's truncation (see the nmod ring);
   flat so that sub-block transforms are simply fewer iterations */
static void
_mtpoly_scale_bias(const sd_fft_ctx_struct * Q, double * d,
                   ulong m_, double b, ulong npts)
{
    vec8d m = vec8d_set_d(vec1d_reduce_0n_to_pmhn((slong) m_, Q->p));
    vec8d bb = vec8d_set_d(b);
    vec8d n = vec8d_set_d(Q->p);
    vec8d ninv = vec8d_set_d(Q->pinv);

    ulong j = 0; do {
        vec8d x0, x1;
        x0 = vec8d_load(d + j + 0);
        x1 = vec8d_load(d + j + 8);
        x0 = vec8d_add(vec8d_mulmod(x0, m, n, ninv), bb);
        x1 = vec8d_add(vec8d_mulmod(x1, m, n, ninv), bb);
        vec8d_store(d + j + 0, x0);
        vec8d_store(d + j + 8, x1);
    } while (j += 16, j < npts);
}

static int
mtpoly_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    mtpoly_ctx_struct * T = MTPOLY_CTX(ctx);
    int status = GR_SUCCESS;
    status |= gr_stream_write(out,
        "Transformed polynomials (fft_small) over multi-word integers mod n"
        " (len_bound ");
    status |= gr_stream_write_si(out, T->N);
    status |= gr_stream_write(out, ", terms_bound ");
    status |= gr_stream_write_ui(out, T->terms_bound);
    status |= gr_stream_write(out, ")");
    return status;
}

static void
mtpoly_ctx_clear(gr_ctx_t ctx)
{
    mtpoly_ctx_struct * T = MTPOLY_CTX(ctx);
    fft_small_plan_clear(T->P);
    flint_free(T->bias_mod_n);
    flint_free(T);
}

static void
mtpoly_init(gr_ptr x, gr_ctx_t ctx)
{
    fft_small_op_init(&MTPOLY(x)->op, MTPOLY_CTX(ctx)->P);
    MTPOLY(x)->len = 0;
    MTPOLY(x)->terms = 0;
    MTPOLY(x)->depth = 0;
    MTPOLY(x)->negs = 0;
}

static void
mtpoly_clear(gr_ptr x, gr_ctx_t FLINT_UNUSED(ctx))
{
    fft_small_op_clear(&MTPOLY(x)->op);
}

static void
mtpoly_swap(gr_ptr x, gr_ptr y, gr_ctx_t FLINT_UNUSED(ctx))
{
    FLINT_SWAP(mtpoly_struct, *MTPOLY(x), *MTPOLY(y));
}

static int
mtpoly_zero(gr_ptr x, gr_ctx_t FLINT_UNUSED(ctx))
{
    MTPOLY(x)->len = 0;
    MTPOLY(x)->terms = 0;
    MTPOLY(x)->depth = 0;
    MTPOLY(x)->negs = 0;
    return GR_SUCCESS;
}

static int
mtpoly_set(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    mtpoly_ctx_struct * T = MTPOLY_CTX(ctx);

    if (res == x)
        return GR_SUCCESS;

    if (MTPOLY(x)->len > 0)
        memcpy(MTPOLY(res)->op.data, MTPOLY(x)->op.data,
               T->P->np * T->P->stride * sizeof(double));
    MTPOLY(res)->op.domain = MTPOLY(x)->op.domain;
    MTPOLY(res)->op.itrunc = MTPOLY(x)->op.itrunc;
    MTPOLY(res)->len = MTPOLY(x)->len;
    MTPOLY(res)->terms = MTPOLY(x)->terms;
    MTPOLY(res)->depth = MTPOLY(x)->depth;
    MTPOLY(res)->negs = MTPOLY(x)->negs;
    return GR_SUCCESS;
}

static truth_t
mtpoly_is_zero(gr_srcptr x, gr_ctx_t FLINT_UNUSED(ctx))
{
    return (MTPOLY(x)->len == 0) ? T_TRUE : T_UNKNOWN;
}

static truth_t
mtpoly_equal(gr_srcptr x, gr_srcptr y, gr_ctx_t FLINT_UNUSED(ctx))
{
    if (MTPOLY(x)->len == 0 && MTPOLY(y)->len == 0)
        return T_TRUE;
    return T_UNKNOWN;
}

static int
_mtpoly_base_ok(gr_ctx_t base_ctx, const mtpoly_ctx_struct * T)
{
    return base_ctx->which_ring == GR_CTX_MPN_MOD &&
           MPN_MOD_CTX_NLIMBS(base_ctx) == T->nlimbs &&
           mpn_cmp(MPN_MOD_CTX_MODULUS(base_ctx),
                   MPN_MOD_CTX_MODULUS(T->base), T->nlimbs) == 0;
}

static int
mtpoly_set_gr_poly(gr_ptr res, gr_srcptr a, slong len,
                   gr_ctx_t base_ctx, gr_ctx_t ctx)
{
    mtpoly_ctx_struct * T = MTPOLY_CTX(ctx);
    nn_srcptr ac = (nn_srcptr) a;

    if (!_mtpoly_base_ok(base_ctx, T))
        return GR_DOMAIN;

    while (len > 0 && flint_mpn_zero_p(ac + (len - 1) * T->nlimbs, T->nlimbs))
        len--;

    if (len > T->N)
        return GR_DOMAIN;

    if (len == 0)
        return mtpoly_zero(res, ctx);

    fft_small_fft_mpn_mod(&MTPOLY(res)->op, ac, len,
                          _op_trunc(len, T->P), T->base, T->P);
    MTPOLY(res)->len = len;
    MTPOLY(res)->terms = 1;
    MTPOLY(res)->depth = 1;
    MTPOLY(res)->negs = 0;
    return GR_SUCCESS;
}

/* shared conversion-out core: reconstruct coefficients [zl, ub) of x,
   ub <= xlen, directly into z (ub - zl coefficients of nlimbs words,
   bias already removed); sdata is per-call temporary space */
static void
_mtpoly_export(mtpoly_ctx_struct * T, gr_srcptr x, slong zl, slong ub,
               double * sdata, nn_ptr z)
{
    const fft_small_plan_struct * P = T->P;
    slong xlen = MTPOLY(x)->len;
    fft_small_op_struct tmp;
    ulong itr = _op_trunc((ulong) xlen, P);
    slong i;

    tmp = MTPOLY(x)->op;
    tmp.data = sdata;

    memcpy(sdata, MTPOLY(x)->op.data, P->np * P->stride * sizeof(double));
    for (i = 0; i < (slong) P->np; i++)
    {
        sd_fft_ctx_struct * Q = P->ffts + P->offset + i;
        double * d = sdata + P->stride * i;
        double b = MTPOLY(x)->negs ? (double) T->bias_res[i] : 0.0;

        sd_ifft_trunc(Q, d, P->depth, itr);
        _mtpoly_scale_bias(Q, d, T->m_orig[i], b, itr);
    }
    fft_small_export_mpn_mod_range(z, &tmp, (ulong) zl, (ulong) ub,
                                   T->base, P);

    if (MTPOLY(x)->negs)
        for (i = 0; i < ub - zl; i++)
            mpn_mod_sub(z + i * T->nlimbs, z + i * T->nlimbs,
                        T->bias_mod_n, T->base);
}

static double *
_mtpoly_scratch_alloc(const mtpoly_ctx_struct * T)
{
    return flint_aligned_alloc(FLINT_FFT_SMALL_ALIGNMENT,
            n_round_up(T->P->np * T->P->stride * sizeof(double), FLINT_FFT_SMALL_ALIGNMENT));
}

/* the GET_GR_POLY_DESTRUCTIVE method: the inverse transforms, scaling
   and bias run on the element's own evaluations, so no scratch is
   allocated and no transform copy is made; the element may only be
   cleared or fully overwritten afterwards */
static int
mtpoly_get_gr_poly_destructive(gr_ptr c, slong * len, gr_ptr x,
                               gr_ctx_t base_ctx, gr_ctx_t ctx)
{
    mtpoly_ctx_struct * T = MTPOLY_CTX(ctx);
    const fft_small_plan_struct * P = T->P;
    nn_ptr cc = (nn_ptr) c;
    slong xlen = MTPOLY(x)->len;
    ulong itr = _op_trunc((ulong) xlen, P);
    slong i;

    if (!_mtpoly_base_ok(base_ctx, T))
        return GR_DOMAIN;

    if (xlen == 0)
    {
        *len = 0;
        return GR_SUCCESS;
    }

    for (i = 0; i < (slong) P->np; i++)
    {
        sd_fft_ctx_struct * Q = P->ffts + P->offset + i;
        double * d = MTPOLY(x)->op.data + P->stride * i;
        double b = MTPOLY(x)->negs ? (double) T->bias_res[i] : 0.0;

        sd_ifft_trunc(Q, d, P->depth, itr);
        _mtpoly_scale_bias(Q, d, T->m_orig[i], b, itr);
    }
    fft_small_export_mpn_mod_range(cc, &MTPOLY(x)->op, 0, (ulong) xlen,
                                   T->base, P);
    if (MTPOLY(x)->negs)
        for (i = 0; i < xlen; i++)
            mpn_mod_sub(cc + i * T->nlimbs, cc + i * T->nlimbs,
                        T->bias_mod_n, T->base);
    MTPOLY(x)->len = 0;

    while (xlen > 0 &&
           flint_mpn_zero_p(cc + (xlen - 1) * T->nlimbs, T->nlimbs))
        xlen--;
    *len = xlen;
    return GR_SUCCESS;
}

static int
mtpoly_get_gr_poly(gr_ptr c, slong * len, gr_srcptr x,
                   gr_ctx_t base_ctx, gr_ctx_t ctx)
{
    mtpoly_ctx_struct * T = MTPOLY_CTX(ctx);
    nn_ptr cc = (nn_ptr) c;
    slong xlen = MTPOLY(x)->len;

    if (!_mtpoly_base_ok(base_ctx, T))
        return GR_DOMAIN;

    if (xlen == 0)
    {
        *len = 0;
        return GR_SUCCESS;
    }

    {
        double * sdata = _mtpoly_scratch_alloc(T);
        _mtpoly_export(T, x, 0, xlen, sdata, cc);
        flint_aligned_free(sdata);
    }

    while (xlen > 0 && flint_mpn_zero_p(cc + (xlen - 1) * T->nlimbs, T->nlimbs))
        xlen--;
    *len = xlen;
    return GR_SUCCESS;
}

static int
mtpoly_get_gr_poly_window(gr_ptr c, gr_srcptr x, slong zl, slong zh,
                          gr_ctx_t base_ctx, gr_ctx_t ctx)
{
    mtpoly_ctx_struct * T = MTPOLY_CTX(ctx);
    nn_ptr cc = (nn_ptr) c;
    slong xlen = MTPOLY(x)->len;
    slong ub;

    if (!_mtpoly_base_ok(base_ctx, T))
        return GR_DOMAIN;
    if (zl < 0 || zh < zl || zh > T->N)
        return GR_DOMAIN;

    ub = FLINT_MIN(zh, xlen);
    if (ub > zl)
    {
        double * sdata = _mtpoly_scratch_alloc(T);
        _mtpoly_export(T, x, zl, ub, sdata, cc);
        flint_aligned_free(sdata);
    }
    else
        ub = zl;
    flint_mpn_zero(cc + (ub - zl) * T->nlimbs, (zh - ub) * T->nlimbs);
    return GR_SUCCESS;
}

static int
mtpoly_one(gr_ptr x, gr_ctx_t ctx)
{
    mtpoly_ctx_struct * T = MTPOLY_CTX(ctx);
    nn_ptr c = flint_calloc(T->nlimbs, sizeof(ulong));
    c[0] = 1;

    fft_small_fft_mpn_mod(&MTPOLY(x)->op, c, 1, _op_trunc(1, T->P), T->base, T->P);
    flint_free(c);
    MTPOLY(x)->len = 1;
    MTPOLY(x)->terms = 1;
    MTPOLY(x)->depth = 1;
    MTPOLY(x)->negs = 0;
    return GR_SUCCESS;
}

static int
mtpoly_neg(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    mtpoly_ctx_struct * T = MTPOLY_CTX(ctx);

    if (MTPOLY(x)->len == 0)
        return mtpoly_zero(res, ctx);

    fft_small_op_neg(&MTPOLY(res)->op, &MTPOLY(x)->op, T->P);
    MTPOLY(res)->op.domain = FFT_SMALL_OP_PRIMAL;
    MTPOLY(res)->len = MTPOLY(x)->len;
    MTPOLY(res)->terms = MTPOLY(x)->terms;
    MTPOLY(res)->depth = MTPOLY(x)->depth;
    MTPOLY(res)->negs = 1;
    return GR_SUCCESS;
}

static int
mtpoly_add(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
{
    mtpoly_ctx_struct * T = MTPOLY_CTX(ctx);
    ulong terms;

    if (MTPOLY(x)->len == 0)
        return mtpoly_set(res, y, ctx);
    if (MTPOLY(y)->len == 0)
        return mtpoly_set(res, x, ctx);

    terms = MTPOLY(x)->terms + MTPOLY(y)->terms;
    if (terms < MTPOLY(x)->terms || terms > T->terms_bound)
        return GR_UNABLE;

    MTPOLY(x)->op.domain = MTPOLY(y)->op.domain = FFT_SMALL_OP_PRIMAL;
    fft_small_op_add(&MTPOLY(res)->op, &MTPOLY(x)->op, &MTPOLY(y)->op, T->P);
    MTPOLY(res)->op.domain = FFT_SMALL_OP_PRIMAL;
    MTPOLY(res)->len = FLINT_MAX(MTPOLY(x)->len, MTPOLY(y)->len);
    MTPOLY(res)->terms = terms;
    MTPOLY(res)->depth = FLINT_MAX(MTPOLY(x)->depth, MTPOLY(y)->depth);
    MTPOLY(res)->negs = MTPOLY(x)->negs | MTPOLY(y)->negs;
    return GR_SUCCESS;
}

static int
mtpoly_sub(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
{
    mtpoly_ctx_struct * T = MTPOLY_CTX(ctx);
    ulong terms;

    if (MTPOLY(y)->len == 0)
        return mtpoly_set(res, x, ctx);
    if (MTPOLY(x)->len == 0)
        return mtpoly_neg(res, y, ctx);

    terms = MTPOLY(x)->terms + MTPOLY(y)->terms;
    if (terms < MTPOLY(x)->terms || terms > T->terms_bound)
        return GR_UNABLE;

    MTPOLY(x)->op.domain = MTPOLY(y)->op.domain = FFT_SMALL_OP_PRIMAL;
    fft_small_op_sub(&MTPOLY(res)->op, &MTPOLY(x)->op, &MTPOLY(y)->op, T->P);
    MTPOLY(res)->op.domain = FFT_SMALL_OP_PRIMAL;
    MTPOLY(res)->len = FLINT_MAX(MTPOLY(x)->len, MTPOLY(y)->len);
    MTPOLY(res)->terms = terms;
    MTPOLY(res)->depth = FLINT_MAX(MTPOLY(x)->depth, MTPOLY(y)->depth);
    MTPOLY(res)->negs = 1;
    return GR_SUCCESS;
}

static int
_mtpoly_mul_terms(ulong * terms, const mtpoly_struct * x,
                  const mtpoly_struct * y, ulong terms_bound)
{
    ulong hi, t;

    umul_ppmm(hi, t, x->terms, y->terms);
    if (hi != 0)
        return 0;
    umul_ppmm(hi, t, t, (ulong) FLINT_MIN(x->len, y->len));
    if (hi != 0 || t > terms_bound)
        return 0;
    *terms = t;
    return 1;
}

static int
mtpoly_mul(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
{
    mtpoly_ctx_struct * T = MTPOLY_CTX(ctx);
    slong len;
    ulong terms;

    if (MTPOLY(x)->len == 0 || MTPOLY(y)->len == 0)
        return mtpoly_zero(res, ctx);

    len = MTPOLY(x)->len + MTPOLY(y)->len - 1;
    if (len > T->N ||
        !_mtpoly_mul_terms(&terms, MTPOLY(x), MTPOLY(y), T->terms_bound))
        return GR_UNABLE;
    if (MTPOLY(x)->depth + MTPOLY(y)->depth > T->max_depth)
        return GR_UNABLE;

    MTPOLY(x)->op.domain = MTPOLY(y)->op.domain = FFT_SMALL_OP_PRIMAL;
    fft_small_op_mul(&MTPOLY(res)->op, &MTPOLY(x)->op, &MTPOLY(y)->op, T->P);
    MTPOLY(res)->op.domain = FFT_SMALL_OP_PRIMAL;
    MTPOLY(res)->len = len;
    MTPOLY(res)->terms = terms;
    MTPOLY(res)->depth = MTPOLY(x)->depth + MTPOLY(y)->depth;
    MTPOLY(res)->negs = MTPOLY(x)->negs | MTPOLY(y)->negs;
    return GR_SUCCESS;
}

static int
mtpoly_addsubmul(gr_ptr res, gr_srcptr x, gr_srcptr y, int sub, gr_ctx_t ctx)
{
    mtpoly_ctx_struct * T = MTPOLY_CTX(ctx);
    slong len;
    ulong terms, t2;

    if (MTPOLY(x)->len == 0 || MTPOLY(y)->len == 0)
        return GR_SUCCESS;

    /* with res aliasing an operand, the per-call domain bookkeeping
       below cannot mark the same struct as both input and accumulator;
       route through the ring's own guarded multiply and add */
    if (res == x || res == y)
    {
        gr_ptr t;
        int status;
        GR_TMP_INIT(t, ctx);
        status = mtpoly_mul(t, x, y, ctx);
        if (status == GR_SUCCESS)
            status = sub ? mtpoly_sub(res, res, t, ctx)
                         : mtpoly_add(res, res, t, ctx);
        GR_TMP_CLEAR(t, ctx);
        return status;
    }

    if (MTPOLY(res)->len == 0)
    {
        int status = mtpoly_mul(res, x, y, ctx);
        if (status == GR_SUCCESS && sub)
            status = mtpoly_neg(res, res, ctx);
        return status;
    }

    len = MTPOLY(x)->len + MTPOLY(y)->len - 1;
    if (len > T->N ||
        !_mtpoly_mul_terms(&t2, MTPOLY(x), MTPOLY(y), T->terms_bound))
        return GR_UNABLE;
    terms = MTPOLY(res)->terms + t2;
    if (terms < t2 || terms > T->terms_bound)
        return GR_UNABLE;
    if (MTPOLY(x)->depth + MTPOLY(y)->depth > T->max_depth)
        return GR_UNABLE;

    MTPOLY(x)->op.domain = MTPOLY(y)->op.domain = FFT_SMALL_OP_PRIMAL;
    MTPOLY(res)->op.domain = FFT_SMALL_OP_PRODUCT;
    if (sub)
        fft_small_op_submul(&MTPOLY(res)->op, &MTPOLY(x)->op,
                            &MTPOLY(y)->op, T->P);
    else
        fft_small_op_addmul(&MTPOLY(res)->op, &MTPOLY(x)->op,
                            &MTPOLY(y)->op, T->P);
    MTPOLY(res)->op.domain = FFT_SMALL_OP_PRIMAL;
    MTPOLY(res)->len = FLINT_MAX(MTPOLY(res)->len, len);
    MTPOLY(res)->terms = terms;
    MTPOLY(res)->depth = FLINT_MAX(MTPOLY(res)->depth,
                            MTPOLY(x)->depth + MTPOLY(y)->depth);
    MTPOLY(res)->negs = sub ? 1 :
        (MTPOLY(res)->negs | MTPOLY(x)->negs | MTPOLY(y)->negs);
    return GR_SUCCESS;
}

static int
mtpoly_addmul(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
{
    return mtpoly_addsubmul(res, x, y, 0, ctx);
}

static int
mtpoly_submul(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
{
    return mtpoly_addsubmul(res, x, y, 1, ctx);
}

static gr_funcptr __mtpoly_methods[GR_METHOD_TAB_SIZE];
static int __mtpoly_methods_initialized = 0;

/* random depth-1 elements and a value writer, via the stored base
   context: gr_test_ring compliance at conversion cost only */
static int
mtpoly_randtest(gr_ptr res, flint_rand_t state, gr_ctx_t ctx)
{
    mtpoly_ctx_struct * T = MTPOLY_CTX(ctx);
    slong nl = T->nlimbs;
    slong i, j, len = n_randint(state, T->N + 1);
    nn_ptr c;
    nn_srcptr m = MPN_MOD_CTX_MODULUS(T->base);
    int status;

    if (len == 0 || n_randint(state, 8) == 0)
        return gr_zero(res, ctx);

    c = flint_malloc(len * nl * sizeof(ulong));
    for (i = 0; i < len; i++)
    {
        nn_ptr ci = c + i * nl;
        for (j = 0; j < nl; j++)
            ci[j] = n_randtest(state);
        /* keep the value below the modulus: top limb strictly below
           the modulus's top limb */
        ci[nl - 1] = (m[nl - 1] > 0) ? n_randint(state, m[nl - 1]) : 0;
    }
    status = mtpoly_set_gr_poly(res, c, len, T->base, ctx);
    flint_free(c);
    return status;
}

static int
mtpoly_write(gr_stream_t out, gr_srcptr x, gr_ctx_t ctx)
{
    mtpoly_ctx_struct * T = MTPOLY_CTX(ctx);
    slong nl = T->nlimbs;
    nn_ptr c;
    slong i, len = 0;
    int status;
    fmpz_t v;

    /* every element's polynomial length is bounded by the ring's N
       (products beyond it are declined at the operation) */
    c = flint_malloc(FLINT_MAX(T->N, 1) * nl * sizeof(ulong));
    status = mtpoly_get_gr_poly(c, &len, x, T->base, ctx);
    if (status == GR_SUCCESS)
    {
        fmpz_init(v);
        status |= gr_stream_write(out, "[");
        for (i = 0; i < len; i++)
        {
            if (i > 0)
                status |= gr_stream_write(out, ", ");
            fmpz_set_ui_array(v, c + i * nl, nl);
            status |= gr_stream_write_fmpz(out, v);
        }
        status |= gr_stream_write(out, "]");
        fmpz_clear(v);
    }
    flint_free(c);
    return status;
}

static gr_method_tab_input __mtpoly_methods_input[] =
{
    {GR_METHOD_CTX_WRITE,       (gr_funcptr) (void (*)(void)) mtpoly_ctx_write},
    {GR_METHOD_CTX_CLEAR,       (gr_funcptr) (void (*)(void)) mtpoly_ctx_clear},
    {GR_METHOD_INIT,            (gr_funcptr) (void (*)(void)) mtpoly_init},
    {GR_METHOD_CLEAR,           (gr_funcptr) (void (*)(void)) mtpoly_clear},
    {GR_METHOD_SWAP,            (gr_funcptr) (void (*)(void)) mtpoly_swap},
    {GR_METHOD_SET,             (gr_funcptr) (void (*)(void)) mtpoly_set},
    {GR_METHOD_ZERO,            (gr_funcptr) (void (*)(void)) mtpoly_zero},
    {GR_METHOD_ONE,             (gr_funcptr) (void (*)(void)) mtpoly_one},
    {GR_METHOD_WRITE,           (gr_funcptr) (void (*)(void)) mtpoly_write},
    {GR_METHOD_RANDTEST,        (gr_funcptr) (void (*)(void)) mtpoly_randtest},
    {GR_METHOD_IS_ZERO,         (gr_funcptr) (void (*)(void)) mtpoly_is_zero},
    {GR_METHOD_EQUAL,           (gr_funcptr) (void (*)(void)) mtpoly_equal},
    {GR_METHOD_NEG,             (gr_funcptr) (void (*)(void)) mtpoly_neg},
    {GR_METHOD_ADD,             (gr_funcptr) (void (*)(void)) mtpoly_add},
    {GR_METHOD_SUB,             (gr_funcptr) (void (*)(void)) mtpoly_sub},
    {GR_METHOD_MUL,             (gr_funcptr) (void (*)(void)) mtpoly_mul},
    {GR_METHOD_ADDMUL,          (gr_funcptr) (void (*)(void)) mtpoly_addmul},
    {GR_METHOD_SUBMUL,          (gr_funcptr) (void (*)(void)) mtpoly_submul},
    {GR_METHOD_SET_GR_POLY,     (gr_funcptr) (void (*)(void)) mtpoly_set_gr_poly},
    {GR_METHOD_GET_GR_POLY,     (gr_funcptr) (void (*)(void)) mtpoly_get_gr_poly},
    {GR_METHOD_GET_GR_POLY_DESTRUCTIVE,
                                (gr_funcptr) (void (*)(void)) mtpoly_get_gr_poly_destructive},
    {GR_METHOD_GET_GR_POLY_WINDOW, (gr_funcptr) (void (*)(void)) mtpoly_get_gr_poly_window},
    {0,                         (gr_funcptr) (void (*)(void)) NULL},
};

int
_gr_mpn_mod_ctx_init_transformed_poly_repr(gr_ctx_t ctx, gr_ctx_t base,
        slong len_bound, slong terms_bound,
        const struct gr_transformed_poly_workload_struct * workload)
{
    mtpoly_ctx_struct * T;
    slong N = len_bound;
    slong nlimbs;
    ulong nbits, i;

    if (base->which_ring != GR_CTX_MPN_MOD || N < 1 || terms_bound < 1 ||
            (ulong) terms_bound > UWORD_MAX / 8)
        return GR_UNABLE;

    nlimbs = MPN_MOD_CTX_NLIMBS(base);
    nbits = MPN_MOD_CTX_MODULUS_BITS(base);

    T = FLINT_ARRAY_ALLOC(1, mtpoly_ctx_struct);
    T->base = base;
    T->nlimbs = nlimbs;
    T->N = N;
    T->terms_bound = (ulong) terms_bound;
    T->max_depth = 2;

    /* provision the chinese remaindering for products of multiplicative
       depth up to 2 with up to terms_bound elementary products, times
       four: one bit for the signed-value bias and one because the
       reconstruction requires values below half the prime product */
    T->P->R = get_default_mpn_ctx();
    T->P->sign = 0;
    _fft_small_plan_set_window(T->P, 0, (ulong) N, (ulong) N,
                               n_round_up((ulong) N, BLK_SZ), LG_BLK_SZ);
    if (!_fft_small_plan_set_bound(T->P, 4 * (ulong) terms_bound,
                                   2 * nbits, MPN_CTX_NCRTS))
    {
        flint_free(T);
        return GR_UNABLE;
    }
    _fft_small_plan_set_normalizers(T->P);

    /* workload cost model, as in the nmod ring; there is no repacked or
       single-prime fused alternative for multi-word moduli, so the fused
       side always uses the same transform sizes with its own (usually
       equal) prime count */
    {
        const gr_transformed_poly_workload_struct def = { 2, 1, 1, 0, 0, 0 };
        const gr_transformed_poly_workload_struct * wl =
            workload ? workload : &def;

    /* forced initialization (tests): only implementation
       bounds below decline; the profitability model and the
       storage budget are policy and are skipped */
    if (!wl->force)
    {
            double L = (double) _len_trunc((ulong) N);
            double lg = (double) FLINT_BIT_COUNT((ulong) L);
            double np = (double) T->P->np;
            double ni = (double) wl->num_inputs;
            double nm = (double) wl->num_muls;
            double no = (double) wl->num_outputs;
            double npf, rep, fuse, bytes, limit;

            /* per-coefficient conversions cost ~nlimbs here */
            rep = np * L * ((ni + no) * (lg + 2.0 * nlimbs) + nm * 2.0)
                    + no * np * L * 3.0 + 5e4;
            npf = (double) ((2 * nbits + FLINT_BIT_COUNT((ulong) L) + 49) / 50);
            npf = FLINT_MAX(npf, 1.0);
            fuse = nm * (npf * L * (3.0 * lg + 4.0 * nlimbs + 2.0)
                         + npf * L * 3.0);

            /* live elements as declared by the driver; a zero declaration
               derives a conservative count from the workload shape */
            {
                double nlive = wl->num_live > 0 ? (double) wl->num_live
                                                : (ni + no + 2.0);
                bytes = nlive * np *
                        (double) sd_fft_ctx_data_size(T->P->depth) * 8.0;
            }
            limit = wl->mem_limit > 0 ? (double) wl->mem_limit
                        : (double) flint_fft_small_max_transformed_ring_size;

            if (rep >= 0.865 * fuse || bytes > limit)
            {
                fft_small_plan_clear(T->P);
                flint_free(T);
                return GR_UNABLE;
            }
    }
    }

    for (i = 0; i < T->P->np; i++)
    {
        T->m_orig[i] = T->P->m[i];
        T->P->m[i] = 1;
    }

    /* bias C = terms_bound * (n-1)^2 (see the nmod ring) */
    {
        fmpz_t C, nz;
        T->bias_mod_n = FLINT_ARRAY_ALLOC(nlimbs, ulong);

        fmpz_init(nz);
        fmpz_set_ui_array(nz, MPN_MOD_CTX_MODULUS(base), nlimbs);
        fmpz_init(C);
        fmpz_sub_ui(C, nz, 1);
        fmpz_mul(C, C, C);
        fmpz_mul_ui(C, C, (ulong) terms_bound);
        for (i = 0; i < T->P->np; i++)
        {
            sd_fft_ctx_struct * Q = T->P->ffts + T->P->offset + i;
            ulong cop = T->P->np == 1 ? 1 :
                *crt_data_co_prime_red(T->P->crts + T->P->np - 1,
                                       T->P->offset + i);
            ulong Ci = fmpz_fdiv_ui(C, Q->mod.n);
            T->bias_res[i] = nmod_mul(Ci, nmod_inv(cop % Q->mod.n, Q->mod),
                                      Q->mod);
        }
        fmpz_mod(C, C, nz);
        fmpz_get_ui_array(T->bias_mod_n, nlimbs, C);
        fmpz_clear(C);
        fmpz_clear(nz);
    }

    ctx->which_ring = GR_CTX_GR_TRANSFORMED_POLY;
    ctx->sizeof_elem = sizeof(mtpoly_struct);
    ctx->size_limit = WORD_MAX;
    GR_CTX_DATA_AS_PTR(ctx) = T;
    ctx->methods = __mtpoly_methods;

    /* concurrent initializations write identical data, which is the
       accepted pattern for gr method tables */
    if (!__mtpoly_methods_initialized)
    {
        gr_method_tab_init(__mtpoly_methods, __mtpoly_methods_input);
        __mtpoly_methods_initialized = 1;
    }

    return GR_SUCCESS;
}

#else /* FLINT_HAVE_FFT_SMALL */

int
_gr_mpn_mod_ctx_init_transformed_poly_repr(gr_ctx_t FLINT_UNUSED(ctx), gr_ctx_t FLINT_UNUSED(base),
        slong FLINT_UNUSED(len_bound), slong FLINT_UNUSED(terms_bound),
        const struct gr_transformed_poly_workload_struct * FLINT_UNUSED(workload))
{
    return GR_UNABLE;
}

#endif
