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
#include "nmod.h"
#include "nmod_vec.h"
#include "gr.h"
#include "gr_poly.h"

#if FLINT_HAVE_FFT_SMALL

#include "machine_vectors.h"
#include "fft_small.h"

/*
    Ring of nmod polynomials of bounded length in fft_small transformed
    representation (prototype).

    A context wraps an fft_small plan of capacity (len_bound, terms_bound):
    elements represent polynomials of length at most len_bound = N over the
    base nmod ring, and the integer coefficients of any represented value,
    viewed as an expression in freshly converted inputs, accumulate at most
    terms_bound elementary products of base elements per coefficient.

    Elements store, per plan prime, the plain pointwise evaluations of the
    represented polynomial on the transform domain: the plan's normalizers
    are overridden to 1 so that the op-layer pointwise kernels compute
    plain products and sums, which makes the representation closed under
    +, -, * (a graded primal/product distinction would otherwise prevent
    multiplying previously multiplied elements). The saved normalizers
    m_i = inv(cop_i * 2^depth) mod p_i are applied once at conversion back
    to coefficients, which is exactly the scaling the chinese remaindering
    expects after an unnormalized inverse transform.

    The chinese remaindering's final reduction handles reconstructed
    integers of at most three limbs, so representable expressions are
    restricted to multiplicative depth two in base elements (products of
    previously multiplied elements return GR_UNABLE beyond that; callers
    renormalize through a conversion roundtrip where deeper products are
    needed). Single-prime plans computing directly mod n are exact and
    carry no depth restriction.

    Values represented by expressions involving subtraction can have
    negative true integer coefficients, while the chinese remaindering
    reconstructs representatives in [0, prod of primes). Conversion back
    to coefficients therefore adds a fixed bias C = terms_bound * (n-1)^2,
    at least as large as any representable coefficient magnitude, to each
    inverse-transformed coefficient (as the per-prime residue of C in the
    convention the chinese remaindering expects), and subtracts C mod n
    from each reconstructed coefficient. The plan is
    provisioned with an accumulation bound of 2 * terms_bound so that the
    chinese remaindering's declared coefficient bound is exactly 2C, which
    both selects enough primes and calibrates the reduction branch
    preconditions for the biased values.

    Each element carries a virtual length and a terms counter; operations
    return GR_UNABLE when a result would exceed the context capacity
    (wraparound past the transform size, or coefficients past the prime
    product), which composes safely through arbitrary expressions.

    Limitations of the prototype, by design:
    - equal / is_zero return T_UNKNOWN for nonzero-length operands (the
      redundant floating point residues are not canonicalized);
    - a context may be shared by any number of threads: conversions use
      plain per-call temporary space, and the context itself is read-only
      after construction (distinct elements may be operated on
      concurrently; a single element still belongs to one thread at a
      time, as usual);
    - the exported window is always [0, N).
*/

typedef struct
{
    fft_small_plan_t P;
    nmod_t mod;
    slong N;                        /* length capacity */
    ulong terms_bound;              /* accumulation capacity */
    ulong max_depth;                /* multiplicative depth capacity */
    ulong m_orig[MPN_CTX_NCRTS];    /* saved export normalizers */
    ulong bias_res[MPN_CTX_NCRTS];  /* C * inv(cop_i) mod p_i */
    ulong bias_mod_n;               /* C mod n */
} tpoly_ctx_struct;

typedef struct
{
    fft_small_op_struct op;
    slong len;                      /* virtual length; 0 = zero element */
    ulong terms;                    /* elementary product count bound */
    ulong depth;                    /* multiplicative depth in base elements */
    int negs;                       /* true integer coeffs may be negative */
} tpoly_struct;

/* must match FFT_SMALL_MAX_REPACK in nmod_poly/mullow_fft_small.c: below
   this, the fused drivers pack two coefficients per transform slot */
#define TPOLY_MAX_REPACK 23

#define TPOLY_CTX(ctx) ((tpoly_ctx_struct *) GR_CTX_DATA_AS_PTR(ctx))
#define TPOLY(x) ((tpoly_struct *) (x))

/* d = m*d + b over the first nblk blocks: applies the saved normalizer
   (and, for values that may be negative, the reconstruction bias) to the
   inverse-transformed coefficients; scaling commutes with the linear
   inverse transform, and applying it afterwards touches only the
   element's virtual length instead of the full transform */
static void
_tpoly_scale_bias(const sd_fft_ctx_struct * Q, double * d,
                  ulong m_, double b, ulong nblk)
{
    vec8d m = vec8d_set_d(vec1d_reduce_0n_to_pmhn((slong) m_, Q->p));
    vec8d bb = vec8d_set_d(b);
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
            x0 = vec8d_add(vec8d_mulmod(x0, m, n, ninv), bb);
            x1 = vec8d_add(vec8d_mulmod(x1, m, n, ninv), bb);
            vec8d_store(dx + j + 0, x0);
            vec8d_store(dx + j + 8, x1);
        } while (j += 16, j < BLK_SZ);
    }
}

static int
tpoly_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    tpoly_ctx_struct * T = TPOLY_CTX(ctx);
    int status = GR_SUCCESS;
    status |= gr_stream_write(out, "Transformed polynomials (fft_small) over integers mod ");
    status |= gr_stream_write_ui(out, T->mod.n);
    status |= gr_stream_write(out, " (len_bound ");
    status |= gr_stream_write_si(out, T->N);
    status |= gr_stream_write(out, ", terms_bound ");
    status |= gr_stream_write_ui(out, T->terms_bound);
    status |= gr_stream_write(out, ")");
    return status;
}

static void
tpoly_ctx_clear(gr_ctx_t ctx)
{
    tpoly_ctx_struct * T = TPOLY_CTX(ctx);
    fft_small_plan_clear(T->P);
    flint_free(T);
}

static void
tpoly_init(gr_ptr x, gr_ctx_t ctx)
{
    fft_small_op_init(&TPOLY(x)->op, TPOLY_CTX(ctx)->P);
    TPOLY(x)->len = 0;
    TPOLY(x)->terms = 0;
    TPOLY(x)->depth = 0;
    TPOLY(x)->negs = 0;
}

static void
tpoly_clear(gr_ptr x, gr_ctx_t FLINT_UNUSED(ctx))
{
    fft_small_op_clear(&TPOLY(x)->op);
}

static void
tpoly_swap(gr_ptr x, gr_ptr y, gr_ctx_t FLINT_UNUSED(ctx))
{
    FLINT_SWAP(tpoly_struct, *TPOLY(x), *TPOLY(y));
}

static int
tpoly_zero(gr_ptr x, gr_ctx_t FLINT_UNUSED(ctx))
{
    TPOLY(x)->len = 0;
    TPOLY(x)->terms = 0;
    TPOLY(x)->depth = 0;
    TPOLY(x)->negs = 0;
    return GR_SUCCESS;
}

static int
tpoly_set(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    tpoly_ctx_struct * T = TPOLY_CTX(ctx);

    if (res == x)
        return GR_SUCCESS;

    if (TPOLY(x)->len > 0)
        memcpy(TPOLY(res)->op.data, TPOLY(x)->op.data,
                     T->P->np * T->P->stride * sizeof(double));
    TPOLY(res)->op.domain = TPOLY(x)->op.domain;
    TPOLY(res)->op.itrunc = TPOLY(x)->op.itrunc;
    TPOLY(res)->len = TPOLY(x)->len;
    TPOLY(res)->terms = TPOLY(x)->terms;
    TPOLY(res)->depth = TPOLY(x)->depth;
    TPOLY(res)->negs = TPOLY(x)->negs;
    return GR_SUCCESS;
}

static truth_t
tpoly_is_zero(gr_srcptr x, gr_ctx_t FLINT_UNUSED(ctx))
{
    return (TPOLY(x)->len == 0) ? T_TRUE : T_UNKNOWN;
}

static truth_t
tpoly_equal(gr_srcptr x, gr_srcptr y, gr_ctx_t FLINT_UNUSED(ctx))
{
    if (TPOLY(x)->len == 0 && TPOLY(y)->len == 0)
        return T_TRUE;
    return T_UNKNOWN;
}

/* conversion from coefficients: the SET_GR_POLY method */
static int
tpoly_set_gr_poly(gr_ptr res, gr_srcptr a, slong len,
                  gr_ctx_t base_ctx, gr_ctx_t ctx)
{
    tpoly_ctx_struct * T = TPOLY_CTX(ctx);
    nn_srcptr ac = (nn_srcptr) a;

    if (base_ctx->which_ring != GR_CTX_NMOD ||
            NMOD_CTX(base_ctx).n != T->mod.n)
        return GR_DOMAIN;

    while (len > 0 && ac[len - 1] == 0)
        len--;

    if (len > T->N)
        return GR_DOMAIN;

    if (len == 0)
        return tpoly_zero(res, ctx);

    fft_small_fft_nmod(&TPOLY(res)->op, ac, len,
                       n_round_up(len, BLK_SZ), T->mod, T->P);
    TPOLY(res)->len = len;
    TPOLY(res)->terms = 1;
    TPOLY(res)->depth = 1;
    TPOLY(res)->negs = 0;
    return GR_SUCCESS;
}

/* shared conversion-out core: reconstruct coefficients [zl, ub) of x,
   ub <= xlen, directly into z (z receives ub - zl values, bias already
   removed); sdata is per-call temporary space of np * stride doubles */
static void
_tpoly_export(tpoly_ctx_struct * T, gr_srcptr x, slong zl, slong ub,
              double * sdata, nn_ptr z)
{
    const fft_small_plan_struct * P = T->P;
    slong xlen = TPOLY(x)->len;
    fft_small_op_struct tmp;
    ulong itr = n_round_up((ulong) xlen, BLK_SZ);
    slong i;

    tmp = TPOLY(x)->op;
    tmp.data = sdata;

    /* invert up to the element's virtual length (truncating below it
       would violate the inverse transform's zero-tail assumption even
       for a smaller window), apply the saved normalizers and, if the
       values may be negative, the reconstruction bias, then remainder
       the window */
    memcpy(sdata, TPOLY(x)->op.data, P->np * P->stride * sizeof(double));
    for (i = 0; i < (slong) P->np; i++)
    {
        sd_fft_ctx_struct * Q = P->ffts + P->offset + i;
        double * d = sdata + P->stride * i;
        double b = TPOLY(x)->negs ? (double) T->bias_res[i] : 0.0;

        sd_ifft_trunc(Q, d, P->depth, itr);
        _tpoly_scale_bias(Q, d, T->m_orig[i], b, itr / BLK_SZ);
    }
    fft_small_export_nmod_range(z, &tmp, (ulong) zl, (ulong) ub, T->mod, P);

    if (TPOLY(x)->negs)
        for (i = 0; i < ub - zl; i++)
            z[i] = nmod_sub(z[i], T->bias_mod_n, T->mod);
}

/* Destructive counterpart of the windowed conversion: the inverse
   transforms, scaling and bias run directly on the element's own
   evaluations, so no scratch is allocated and no transform-sized copy is
   made. The element is consumed -- afterwards it may only be cleared or
   fully overwritten as the destination of a set or a multiplication.
   Safe under the drivers' threading, since each thread owns the
   accumulator it converts. */
/* invert, scale and bias the element's own evaluations in place */
static void
_tpoly_invert_in_place(tpoly_ctx_struct * T, gr_ptr x)
{
    const fft_small_plan_struct * P = T->P;
    ulong itr = n_round_up((ulong) TPOLY(x)->len, BLK_SZ);
    slong i;

    for (i = 0; i < (slong) P->np; i++)
    {
        sd_fft_ctx_struct * Q = P->ffts + P->offset + i;
        double * d = TPOLY(x)->op.data + P->stride * i;
        double b = TPOLY(x)->negs ? (double) T->bias_res[i] : 0.0;

        sd_ifft_trunc(Q, d, P->depth, itr);
        _tpoly_scale_bias(Q, d, T->m_orig[i], b, itr / BLK_SZ);
    }
}

int
_gr_nmod_tpoly_get_gr_poly_window_destructive(nn_ptr cc, gr_ptr x,
        slong zl, slong zh, gr_ctx_t ctx)
{
    tpoly_ctx_struct * T = TPOLY_CTX(ctx);
    slong xlen = TPOLY(x)->len;
    slong ub, i;

    if (zl < 0 || zh < zl || zh > T->N)
        return GR_DOMAIN;

    ub = FLINT_MIN(zh, xlen);
    if (ub > zl)
    {
        _tpoly_invert_in_place(T, x);
        fft_small_export_nmod_range(cc, &TPOLY(x)->op, (ulong) zl,
                                    (ulong) ub, T->mod, T->P);

        if (TPOLY(x)->negs)
            for (i = 0; i < ub - zl; i++)
                cc[i] = nmod_sub(cc[i], T->bias_mod_n, T->mod);
    }
    else
        ub = zl;
    memset(cc + (ub - zl), 0, (zh - ub) * sizeof(ulong));

    TPOLY(x)->len = 0;
    return GR_SUCCESS;
}

/* the GET_GR_POLY_DESTRUCTIVE method: full-length conversion out on the
   element's own storage */
static int
tpoly_get_gr_poly_destructive(gr_ptr c, slong * len, gr_ptr x,
                              gr_ctx_t base_ctx, gr_ctx_t ctx)
{
    tpoly_ctx_struct * T = TPOLY_CTX(ctx);
    nn_ptr cc = (nn_ptr) c;
    slong xlen = TPOLY(x)->len;
    slong i;

    if (base_ctx->which_ring != GR_CTX_NMOD ||
            NMOD_CTX(base_ctx).n != T->mod.n)
        return GR_DOMAIN;

    if (xlen == 0)
    {
        *len = 0;
        return GR_SUCCESS;
    }

    _tpoly_invert_in_place(T, x);
    fft_small_export_nmod_range(cc, &TPOLY(x)->op, 0, (ulong) xlen,
                                T->mod, T->P);
    if (TPOLY(x)->negs)
        for (i = 0; i < xlen; i++)
            cc[i] = nmod_sub(cc[i], T->bias_mod_n, T->mod);
    TPOLY(x)->len = 0;

    while (xlen > 0 && cc[xlen - 1] == 0)
        xlen--;
    *len = xlen;
    return GR_SUCCESS;
}

/* per-call temporary for the inverse-transform copy */
static double *
_tpoly_scratch_alloc(const tpoly_ctx_struct * T)
{
    return flint_aligned_alloc(FLINT_FFT_SMALL_ALIGNMENT,
            n_round_up(T->P->np * T->P->stride * sizeof(double), FLINT_FFT_SMALL_ALIGNMENT));
}

/* conversion to coefficients: the GET_GR_POLY method; c must have room
   for the element's virtual length; the returned length is normalized */
static int
tpoly_get_gr_poly(gr_ptr c, slong * len, gr_srcptr x,
                  gr_ctx_t base_ctx, gr_ctx_t ctx)
{
    tpoly_ctx_struct * T = TPOLY_CTX(ctx);
    nn_ptr cc = (nn_ptr) c;
    slong xlen = TPOLY(x)->len;

    if (base_ctx->which_ring != GR_CTX_NMOD ||
            NMOD_CTX(base_ctx).n != T->mod.n)
        return GR_DOMAIN;

    if (xlen == 0)
    {
        *len = 0;
        return GR_SUCCESS;
    }

    {
        double * sdata = _tpoly_scratch_alloc(T);
        _tpoly_export(T, x, 0, xlen, sdata, cc);
        flint_aligned_free(sdata);
    }

    while (xlen > 0 && cc[xlen - 1] == 0)
        xlen--;
    *len = xlen;
    return GR_SUCCESS;
}

/* windowed conversion: writes the coefficients [zl, zh) of the
   represented polynomial into c (zeros beyond its virtual length) */
static int
tpoly_get_gr_poly_window(gr_ptr c, gr_srcptr x, slong zl, slong zh,
                         gr_ctx_t base_ctx, gr_ctx_t ctx)
{
    tpoly_ctx_struct * T = TPOLY_CTX(ctx);
    nn_ptr cc = (nn_ptr) c;
    slong xlen = TPOLY(x)->len;
    slong ub;

    if (base_ctx->which_ring != GR_CTX_NMOD ||
            NMOD_CTX(base_ctx).n != T->mod.n)
        return GR_DOMAIN;
    if (zl < 0 || zh < zl || zh > T->N)
        return GR_DOMAIN;

    ub = FLINT_MIN(zh, xlen);
    if (ub > zl)
    {
        double * sdata = _tpoly_scratch_alloc(T);
        _tpoly_export(T, x, zl, ub, sdata, cc);
        flint_aligned_free(sdata);
    }
    else
        ub = zl;
    memset(cc + (ub - zl), 0, (zh - ub) * sizeof(ulong));
    return GR_SUCCESS;
}

static int
tpoly_one(gr_ptr x, gr_ctx_t ctx)
{
    tpoly_ctx_struct * T = TPOLY_CTX(ctx);
    ulong c = 1;

    fft_small_fft_nmod(&TPOLY(x)->op, &c, 1, BLK_SZ, T->mod, T->P);
    TPOLY(x)->len = 1;
    TPOLY(x)->terms = 1;
    TPOLY(x)->depth = 1;
    TPOLY(x)->negs = 0;
    return GR_SUCCESS;
}

static int
tpoly_neg(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    tpoly_ctx_struct * T = TPOLY_CTX(ctx);

    if (TPOLY(x)->len == 0)
        return tpoly_zero(res, ctx);

    fft_small_op_neg(&TPOLY(res)->op, &TPOLY(x)->op, T->P);
    TPOLY(res)->op.domain = FFT_SMALL_OP_PRIMAL;
    TPOLY(res)->len = TPOLY(x)->len;
    TPOLY(res)->terms = TPOLY(x)->terms;
    TPOLY(res)->depth = TPOLY(x)->depth;
    TPOLY(res)->negs = 1;
    return GR_SUCCESS;
}

static int
tpoly_add(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
{
    tpoly_ctx_struct * T = TPOLY_CTX(ctx);
    ulong terms;

    if (TPOLY(x)->len == 0)
        return tpoly_set(res, y, ctx);
    if (TPOLY(y)->len == 0)
        return tpoly_set(res, x, ctx);

    terms = TPOLY(x)->terms + TPOLY(y)->terms;
    if (terms < TPOLY(x)->terms || terms > T->terms_bound)
        return GR_UNABLE;

    TPOLY(x)->op.domain = TPOLY(y)->op.domain = FFT_SMALL_OP_PRIMAL;
    fft_small_op_add(&TPOLY(res)->op, &TPOLY(x)->op, &TPOLY(y)->op, T->P);
    TPOLY(res)->op.domain = FFT_SMALL_OP_PRIMAL;
    TPOLY(res)->len = FLINT_MAX(TPOLY(x)->len, TPOLY(y)->len);
    TPOLY(res)->terms = terms;
    TPOLY(res)->depth = FLINT_MAX(TPOLY(x)->depth, TPOLY(y)->depth);
    TPOLY(res)->negs = TPOLY(x)->negs | TPOLY(y)->negs;
    return GR_SUCCESS;
}

static int
tpoly_sub(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
{
    tpoly_ctx_struct * T = TPOLY_CTX(ctx);
    ulong terms;

    if (TPOLY(y)->len == 0)
        return tpoly_set(res, x, ctx);
    if (TPOLY(x)->len == 0)
        return tpoly_neg(res, y, ctx);

    terms = TPOLY(x)->terms + TPOLY(y)->terms;
    if (terms < TPOLY(x)->terms || terms > T->terms_bound)
        return GR_UNABLE;

    TPOLY(x)->op.domain = TPOLY(y)->op.domain = FFT_SMALL_OP_PRIMAL;
    fft_small_op_sub(&TPOLY(res)->op, &TPOLY(x)->op, &TPOLY(y)->op, T->P);
    TPOLY(res)->op.domain = FFT_SMALL_OP_PRIMAL;
    TPOLY(res)->len = FLINT_MAX(TPOLY(x)->len, TPOLY(y)->len);
    TPOLY(res)->terms = terms;
    TPOLY(res)->depth = FLINT_MAX(TPOLY(x)->depth, TPOLY(y)->depth);
    TPOLY(res)->negs = 1;
    return GR_SUCCESS;
}

/* terms bound of a product: tx * ty * min(lenx, leny), saturating */
static int
_mul_terms(ulong * terms, const tpoly_struct * x, const tpoly_struct * y,
           ulong terms_bound)
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
tpoly_mul(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
{
    tpoly_ctx_struct * T = TPOLY_CTX(ctx);
    slong len;
    ulong terms;

    if (TPOLY(x)->len == 0 || TPOLY(y)->len == 0)
        return tpoly_zero(res, ctx);

    len = TPOLY(x)->len + TPOLY(y)->len - 1;
    if (len > T->N || !_mul_terms(&terms, TPOLY(x), TPOLY(y), T->terms_bound))
        return GR_UNABLE;
    if (TPOLY(x)->depth + TPOLY(y)->depth > T->max_depth)
        return GR_UNABLE;

    TPOLY(x)->op.domain = TPOLY(y)->op.domain = FFT_SMALL_OP_PRIMAL;
    fft_small_op_mul(&TPOLY(res)->op, &TPOLY(x)->op, &TPOLY(y)->op, T->P);
    TPOLY(res)->op.domain = FFT_SMALL_OP_PRIMAL;
    TPOLY(res)->len = len;
    TPOLY(res)->terms = terms;
    TPOLY(res)->depth = TPOLY(x)->depth + TPOLY(y)->depth;
    TPOLY(res)->negs = TPOLY(x)->negs | TPOLY(y)->negs;
    return GR_SUCCESS;
}

static int
tpoly_addsubmul(gr_ptr res, gr_srcptr x, gr_srcptr y, int sub, gr_ctx_t ctx)
{
    tpoly_ctx_struct * T = TPOLY_CTX(ctx);
    slong len;
    ulong terms, t2;

    if (TPOLY(x)->len == 0 || TPOLY(y)->len == 0)
        return GR_SUCCESS;

    /* with res aliasing an operand, the per-call domain bookkeeping
       below cannot mark the same struct as both input and accumulator;
       route through the ring's own guarded multiply and add */
    if (res == x || res == y)
    {
        gr_ptr t;
        int status;
        GR_TMP_INIT(t, ctx);
        status = tpoly_mul(t, x, y, ctx);
        if (status == GR_SUCCESS)
            status = sub ? tpoly_sub(res, res, t, ctx)
                         : tpoly_add(res, res, t, ctx);
        GR_TMP_CLEAR(t, ctx);
        return status;
    }

    if (TPOLY(res)->len == 0)
    {
        int status = tpoly_mul(res, x, y, ctx);
        if (status == GR_SUCCESS && sub)
            status = tpoly_neg(res, res, ctx);
        return status;
    }

    len = TPOLY(x)->len + TPOLY(y)->len - 1;
    if (len > T->N || !_mul_terms(&t2, TPOLY(x), TPOLY(y), T->terms_bound))
        return GR_UNABLE;
    terms = TPOLY(res)->terms + t2;
    if (terms < t2 || terms > T->terms_bound)
        return GR_UNABLE;
    if (TPOLY(x)->depth + TPOLY(y)->depth > T->max_depth)
        return GR_UNABLE;

    TPOLY(x)->op.domain = TPOLY(y)->op.domain = FFT_SMALL_OP_PRIMAL;
    TPOLY(res)->op.domain = FFT_SMALL_OP_PRODUCT;
    if (sub)
        fft_small_op_submul(&TPOLY(res)->op, &TPOLY(x)->op, &TPOLY(y)->op, T->P);
    else
        fft_small_op_addmul(&TPOLY(res)->op, &TPOLY(x)->op, &TPOLY(y)->op, T->P);
    TPOLY(res)->op.domain = FFT_SMALL_OP_PRIMAL;
    TPOLY(res)->len = FLINT_MAX(TPOLY(res)->len, len);
    TPOLY(res)->terms = terms;
    TPOLY(res)->depth = FLINT_MAX(TPOLY(res)->depth,
                            TPOLY(x)->depth + TPOLY(y)->depth);
    TPOLY(res)->negs = sub ? 1 :
        (TPOLY(res)->negs | TPOLY(x)->negs | TPOLY(y)->negs);
    return GR_SUCCESS;
}

static int
tpoly_addmul(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
{
    return tpoly_addsubmul(res, x, y, 0, ctx);
}

static int
tpoly_submul(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
{
    return tpoly_addsubmul(res, x, y, 1, ctx);
}

static gr_funcptr __tpoly_methods[GR_METHOD_TAB_SIZE];
static int __tpoly_methods_initialized = 0;

/* random depth-1 elements within the length capacity, and a value
   writer via the nondestructive conversion out: what gr_test_ring and
   diagnostics need, at conversion cost only they pay */
static int
tpoly_randtest(gr_ptr res, flint_rand_t state, gr_ctx_t ctx)
{
    tpoly_ctx_struct * T = TPOLY_CTX(ctx);
    gr_ctx_t base;
    slong i, len = n_randint(state, T->N + 1);
    nn_ptr c;
    int status;

    if (len == 0 || n_randint(state, 8) == 0)
        return gr_zero(res, ctx);

    gr_ctx_init_nmod(base, T->mod.n);
    c = flint_malloc(len * sizeof(ulong));
    for (i = 0; i < len; i++)
        c[i] = n_randint(state, T->mod.n);
    c[len - 1] += (c[len - 1] == 0 && T->mod.n > 1);
    status = tpoly_set_gr_poly(res, c, len, base, ctx);
    flint_free(c);
    gr_ctx_clear(base);
    return status;
}

static int
tpoly_write(gr_stream_t out, gr_srcptr x, gr_ctx_t ctx)
{
    tpoly_ctx_struct * T = TPOLY_CTX(ctx);
    gr_ctx_t base;
    nn_ptr c;
    slong i, len = 0;
    int status;

    gr_ctx_init_nmod(base, T->mod.n);
    /* every element's polynomial length is bounded by the ring's N
       (products beyond it are declined at the operation) */
    c = flint_malloc(FLINT_MAX(T->N, 1) * sizeof(ulong));
    status = tpoly_get_gr_poly(c, &len, x, base, ctx);
    if (status == GR_SUCCESS)
    {
        status |= gr_stream_write(out, "[");
        for (i = 0; i < len; i++)
        {
            if (i > 0)
                status |= gr_stream_write(out, ", ");
            status |= gr_stream_write_ui(out, c[i]);
        }
        status |= gr_stream_write(out, "]");
    }
    flint_free(c);
    gr_ctx_clear(base);
    return status;
}

static gr_method_tab_input __tpoly_methods_input[] =
{
    {GR_METHOD_CTX_WRITE,       (gr_funcptr) (void (*)(void)) tpoly_ctx_write},
    {GR_METHOD_CTX_CLEAR,       (gr_funcptr) (void (*)(void)) tpoly_ctx_clear},
    {GR_METHOD_INIT,            (gr_funcptr) (void (*)(void)) tpoly_init},
    {GR_METHOD_CLEAR,           (gr_funcptr) (void (*)(void)) tpoly_clear},
    {GR_METHOD_SWAP,            (gr_funcptr) (void (*)(void)) tpoly_swap},
    {GR_METHOD_SET,             (gr_funcptr) (void (*)(void)) tpoly_set},
    {GR_METHOD_ZERO,            (gr_funcptr) (void (*)(void)) tpoly_zero},
    {GR_METHOD_ONE,             (gr_funcptr) (void (*)(void)) tpoly_one},
    {GR_METHOD_WRITE,           (gr_funcptr) (void (*)(void)) tpoly_write},
    {GR_METHOD_RANDTEST,        (gr_funcptr) (void (*)(void)) tpoly_randtest},
    {GR_METHOD_IS_ZERO,         (gr_funcptr) (void (*)(void)) tpoly_is_zero},
    {GR_METHOD_EQUAL,           (gr_funcptr) (void (*)(void)) tpoly_equal},
    {GR_METHOD_NEG,             (gr_funcptr) (void (*)(void)) tpoly_neg},
    {GR_METHOD_ADD,             (gr_funcptr) (void (*)(void)) tpoly_add},
    {GR_METHOD_SUB,             (gr_funcptr) (void (*)(void)) tpoly_sub},
    {GR_METHOD_MUL,             (gr_funcptr) (void (*)(void)) tpoly_mul},
    {GR_METHOD_ADDMUL,          (gr_funcptr) (void (*)(void)) tpoly_addmul},
    {GR_METHOD_SUBMUL,          (gr_funcptr) (void (*)(void)) tpoly_submul},
    {GR_METHOD_SET_GR_POLY,     (gr_funcptr) (void (*)(void)) tpoly_set_gr_poly},
    {GR_METHOD_GET_GR_POLY,     (gr_funcptr) (void (*)(void)) tpoly_get_gr_poly},
    {GR_METHOD_GET_GR_POLY_DESTRUCTIVE,
                                (gr_funcptr) (void (*)(void)) tpoly_get_gr_poly_destructive},
    {GR_METHOD_GET_GR_POLY_WINDOW, (gr_funcptr) (void (*)(void)) tpoly_get_gr_poly_window},
    {0,                         (gr_funcptr) (void (*)(void)) NULL},
};

/* the nmod overload of gr_ctx_init_gr_poly_transformed_repr */
int
_gr_nmod_ctx_init_transformed_poly_repr(gr_ctx_t ctx, gr_ctx_t base,
                                        slong len_bound, slong terms_bound,
                                        const gr_transformed_poly_workload_struct * workload)
{
    tpoly_ctx_struct * T;
    nmod_t mod;
    slong N = len_bound;
    ulong modbits, i;

    if (base->which_ring != GR_CTX_NMOD || N < 1 || terms_bound < 1 ||
            (ulong) terms_bound > UWORD_MAX / 4)
        return GR_UNABLE;

    mod = NMOD_CTX(base);
    modbits = FLINT_BITS - mod.norm;

    T = FLINT_ARRAY_ALLOC(1, tpoly_ctx_struct);
    T->mod = mod;
    T->N = N;
    T->terms_bound = (ulong) terms_bound;

    /* provision the chinese remaindering for products of multiplicative
       depth up to 2 in base elements with up to terms_bound elementary
       products, doubled to cover the signed-value bias; the crt's final
       reduction handles three limbs, which bounds the total */
    if (2 * modbits + 1 + FLINT_BIT_COUNT((ulong) terms_bound) > 3 * FLINT_BITS
        || !fft_small_plan_init_nmod(T->P, get_default_mpn_ctx(),
            0, N, N, n_round_up(N, BLK_SZ), 2 * (ulong) terms_bound,
            2 * modbits, mod, N))
    {
        flint_free(T);
        return GR_UNABLE;
    }

    /* single-prime plans that compute directly mod n (an fft prime as
       the modulus, or a direct transform) are exact for any expression:
       no depth or accumulation restriction applies, and no signed-value
       bias is needed */
    int exact = T->P->use_direct_fft || (T->P->np == 1 &&
            T->P->ffts[T->P->offset].mod.n == mod.n);

    /* Profitability and memory model. Estimate the cost of performing
       the declared workload in this representation against fused
       multiplications, in double-ops per coefficient: a transform costs
       ~lg(L) per coefficient per prime, a conversion pass ~2, a pointwise
       multiplication ~2, and chinese remaindering ~3. Fused products use
       tighter operand lengths in typical algorithmic callers (factor
       ~0.6), often need fewer primes than the level-wide bound here
       requires, and repack two coefficients per slot for tiny moduli.
       Constants are calibrated against gcd measurements on this class of
       hardware; the model errs conservative. */
    {
        const gr_transformed_poly_workload_struct def = { 2, 1, 1, 0, 0, 0 };
        const gr_transformed_poly_workload_struct * wl =
            workload ? workload : &def;

    /* forced initialization (tests): only implementation
       bounds below decline; the profitability model and the
       storage budget are policy and are skipped */
    if (!wl->force)
    {
            double L = (double) n_round_up(N, BLK_SZ);
            double lg = (double) FLINT_BIT_COUNT((ulong) L);
            double np = (double) T->P->np;
            double ni = (double) wl->num_inputs;
            double nm = (double) wl->num_muls;
            double no = (double) wl->num_outputs;
            double npf, Lf, lgf, rep, fuse, margin, bytes, limit;

            rep = np * L * ((ni + no) * (lg + 2.0) + nm * 2.0)
                    + no * np * L * 3.0 + 5e4;

            if (mod.n <= TPOLY_MAX_REPACK)
            {
                npf = 1.0;
                Lf = L / 2;
            }
            else
            {
                npf = (double) ((2 * modbits + FLINT_BIT_COUNT((ulong) L) + 49)
                                / 50);
                npf = FLINT_MAX(npf, 1.0);
                Lf = L;
            }
            lgf = (double) FLINT_BIT_COUNT((ulong) Lf);
            /* homogeneous workloads (batched operands of comparable length,
               as in matrix products) face fused alternatives at the same
               transform sizes; heterogeneous whole-algorithm batching proved
               unprofitable (operands transformed at the shared upper-bound
               length lose the reuse advantage) and is not attempted */
            fuse = nm * (npf * Lf * (3.0 * lgf + 6.0) + npf * Lf * 3.0);

            /* single-prime transforms leave little to deduplicate relative to
               the surrounding fixed and conversion costs (measured), so
               require a decisive advantage there */
            margin = (T->P->np == 1) ? 0.6 : 0.865;

            /* low-reuse workloads whose live elements far exceed cache go
               memory bound at one or two primes (measured); high-reuse
               workloads amortize the streaming */
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

            if (rep >= margin * fuse
                || (T->P->np <= 2 && nm < 4.0 * ni && bytes > 32.0 * 1048576.0)
                || bytes > limit)
            {
                fft_small_plan_clear(T->P);
                flint_free(T);
                return GR_UNABLE;
            }
    }
    }

    if (exact)
    {
        T->max_depth = UWORD_MAX / 4;
        T->terms_bound = UWORD_MAX / 4;
    }
    else
        T->max_depth = 2;

    /* elements hold plain evaluations: make the pointwise kernels
       apply no scaling, and save the true normalizers for conversion
       back to coefficients */
    for (i = 0; i < T->P->np; i++)
    {
        T->m_orig[i] = T->P->m[i];
        T->P->m[i] = 1;
    }


    /* per-prime residues of the bias C = terms_bound * (n-1)^2, in the
       convention of inverse-transformed coefficients (divided by the crt
       cofactor), to be added to each coefficient before reconstruction;
       exact plans need no bias since all arithmetic is mod n */
    if (exact)
    {
        for (i = 0; i < T->P->np; i++)
            T->bias_res[i] = 0;
        T->bias_mod_n = 0;
    }
    else
    {
        fmpz_t C;

        fmpz_init_set_ui(C, mod.n - 1);
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
        T->bias_mod_n = fmpz_fdiv_ui(C, mod.n);
        fmpz_clear(C);
    }

    ctx->which_ring = GR_CTX_GR_TRANSFORMED_POLY;
    ctx->sizeof_elem = sizeof(tpoly_struct);
    ctx->size_limit = WORD_MAX;
    GR_CTX_DATA_AS_PTR(ctx) = T;
    ctx->methods = __tpoly_methods;

    /* concurrent initializations write identical data, which is the
       accepted pattern for gr method tables */
    if (!__tpoly_methods_initialized)
    {
        gr_method_tab_init(__tpoly_methods, __tpoly_methods_input);
        __tpoly_methods_initialized = 1;
    }

    return GR_SUCCESS;
}

#else /* FLINT_HAVE_FFT_SMALL */

int
_gr_nmod_ctx_init_transformed_poly_repr(gr_ctx_t FLINT_UNUSED(ctx),
        gr_ctx_t FLINT_UNUSED(base), slong FLINT_UNUSED(len_bound),
        slong FLINT_UNUSED(terms_bound),
        const gr_transformed_poly_workload_struct * FLINT_UNUSED(workload))
{
    return GR_UNABLE;
}

#endif
