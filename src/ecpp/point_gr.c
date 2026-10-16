/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5.1

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mod.h"
#include "gr.h"
#include "mpn_mod.h"
#include "ecpp.h"

/*
    Point arithmetic on y^2 = x^3 + a x + b over a generic ring (gr): the
    same Jacobian formulas as point.c, so that the modular arithmetic can
    be mpn_mod for moduli of up to MPN_MOD_MAX_LIMBS limbs and fmpz_mod
    beyond. A point is three consecutive ring elements (X, Y, Z); Z = 0 is
    the point at infinity. All operations are in a ring where every
    nonzero element is expected to be a unit (n prime); the quantities
    whose invertibility is assumed are multiplied into acc, and the
    functions return 0 as soon as acc is zero.
*/

#define PX(P) (P)
#define PY(P) GR_ENTRY(P, 1, sz)
#define PZ(P) GR_ENTRY(P, 2, sz)

#define GR_MUST(expr) do { if ((expr) != GR_SUCCESS) return 0; } while (0)

static int
_gpoint_set(gr_ptr R, gr_srcptr P, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    GR_MUST(gr_set(PX(R), PX(P), ctx));
    GR_MUST(gr_set(PY(R), PY(P), ctx));
    GR_MUST(gr_set(PZ(R), PZ(P), ctx));
    return 1;
}

static int
_gpoint_double(gr_ptr R, gr_srcptr P, gr_srcptr a, gr_ptr acc,
        gr_ptr t1, gr_ptr t2, gr_ptr t3, gr_ptr t4, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;

    if (gr_is_zero(PZ(P), ctx) == T_TRUE || gr_is_zero(PY(P), ctx) == T_TRUE)
    {
        GR_MUST(gr_zero(PZ(R), ctx));
        return 1;
    }
    GR_MUST(gr_mul(acc, acc, PY(P), ctx));

    /* t1 = M = 3 X^2 + a Z^4 */
    GR_MUST(gr_sqr(t1, PX(P), ctx));
    GR_MUST(gr_add(t2, t1, t1, ctx));
    GR_MUST(gr_add(t1, t1, t2, ctx));
    GR_MUST(gr_sqr(t2, PZ(P), ctx));
    GR_MUST(gr_sqr(t2, t2, ctx));
    GR_MUST(gr_mul(t2, t2, a, ctx));
    GR_MUST(gr_add(t1, t1, t2, ctx));

    /* t2 = S = 4 X Y^2, t3 = 8 Y^4 */
    GR_MUST(gr_sqr(t3, PY(P), ctx));
    GR_MUST(gr_mul(t2, PX(P), t3, ctx));
    GR_MUST(gr_add(t2, t2, t2, ctx));
    GR_MUST(gr_add(t2, t2, t2, ctx));
    GR_MUST(gr_sqr(t3, t3, ctx));
    GR_MUST(gr_add(t3, t3, t3, ctx));
    GR_MUST(gr_add(t3, t3, t3, ctx));
    GR_MUST(gr_add(t3, t3, t3, ctx));

    /* Z3 = 2 Y Z */
    GR_MUST(gr_mul(t4, PY(P), PZ(P), ctx));
    GR_MUST(gr_add(PZ(R), t4, t4, ctx));

    /* X3 = M^2 - 2S */
    GR_MUST(gr_sqr(t4, t1, ctx));
    GR_MUST(gr_sub(t4, t4, t2, ctx));
    GR_MUST(gr_sub(PX(R), t4, t2, ctx));

    /* Y3 = M (S - X3) - 8 Y^4 */
    GR_MUST(gr_sub(t2, t2, PX(R), ctx));
    GR_MUST(gr_mul(t2, t2, t1, ctx));
    GR_MUST(gr_sub(PY(R), t2, t3, ctx));
    return 1;
}

/* R = P + Q with Q affine (Z = 1) */
static int
_gpoint_add_affine(gr_ptr R, gr_srcptr P, gr_srcptr Q, gr_srcptr a, gr_ptr acc,
        gr_ptr t1, gr_ptr t2, gr_ptr t3, gr_ptr t4, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;

    if (gr_is_zero(PZ(P), ctx) == T_TRUE)
        return _gpoint_set(R, Q, ctx);

    /* t1 = Z1^2, t2 = U2 = X2 Z1^2, t3 = S2 = Y2 Z1^3 */
    GR_MUST(gr_sqr(t1, PZ(P), ctx));
    GR_MUST(gr_mul(t2, PX(Q), t1, ctx));
    GR_MUST(gr_mul(t1, t1, PZ(P), ctx));
    GR_MUST(gr_mul(t3, PY(Q), t1, ctx));

    /* t2 = H = U2 - X1, t3 = r = S2 - Y1 */
    GR_MUST(gr_sub(t2, t2, PX(P), ctx));
    GR_MUST(gr_sub(t3, t3, PY(P), ctx));

    if (gr_is_zero(t2, ctx) == T_TRUE)
    {
        if (gr_is_zero(t3, ctx) == T_TRUE)
            return _gpoint_double(R, P, a, acc, t1, t2, t3, t4, ctx);
        GR_MUST(gr_zero(PZ(R), ctx));   /* P = -Q */
        return 1;
    }

    GR_MUST(gr_mul(acc, acc, t2, ctx));

    /* Z3 = Z1 H */
    GR_MUST(gr_mul(PZ(R), PZ(P), t2, ctx));

    /* t1 = H^2, t4 = H^3, t1 = V = X1 H^2 */
    GR_MUST(gr_sqr(t1, t2, ctx));
    GR_MUST(gr_mul(t4, t1, t2, ctx));
    GR_MUST(gr_mul(t1, t1, PX(P), ctx));

    /* X3 = r^2 - H^3 - 2V */
    GR_MUST(gr_sqr(t2, t3, ctx));
    GR_MUST(gr_sub(t2, t2, t4, ctx));
    GR_MUST(gr_sub(t2, t2, t1, ctx));
    GR_MUST(gr_sub(PX(R), t2, t1, ctx));

    /* Y3 = r (V - X3) - Y1 H^3 */
    GR_MUST(gr_sub(t1, t1, PX(R), ctx));
    GR_MUST(gr_mul(t1, t1, t3, ctx));
    GR_MUST(gr_mul(t4, t4, PY(P), ctx));
    GR_MUST(gr_sub(PY(R), t1, t4, ctx));
    return 1;
}

/* plain double-and-add; P affine */
static int
_gpoint_mul_binary(gr_ptr R, gr_srcptr P, const fmpz_t k, gr_srcptr a, gr_ptr acc,
        gr_ptr t1, gr_ptr t2, gr_ptr t3, gr_ptr t4, gr_ctx_t ctx)
{
    slong i, bits = fmpz_bits(k);
    gr_ptr T = gr_heap_init_vec(3, ctx);
    int ok = 1;

    ok = _gpoint_set(T, P, ctx);
    for (i = bits - 2; i >= 0 && ok; i--)
    {
        ok = _gpoint_double(T, T, a, acc, t1, t2, t3, t4, ctx);
        if (ok && fmpz_tstbit(k, i))
            ok = _gpoint_add_affine(T, T, P, a, acc, t1, t2, t3, t4, ctx);
        if (gr_is_zero(acc, ctx) == T_TRUE)
            ok = 0;
    }
    if (ok)
        ok = _gpoint_set(R, T, ctx);
    gr_heap_clear_vec(T, 3, ctx);
    return ok;
}

/*
    Normalises num Jacobian points to affine coordinates with one
    inversion, multiplying acc by the Z's. Returns 0 if some Z is zero or
    the product is not invertible (the points are then left alone).
*/
static int
_gpoints_to_affine(gr_ptr T, slong num, gr_ptr acc, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem, i;
    gr_ptr prod = gr_heap_init_vec(num, ctx);
    gr_ptr inv, t, g;
    int ok = 1;

    GR_TMP_INIT3(inv, t, g, ctx);
    for (i = 0; i < num && ok; i++)
        if (gr_is_zero(PZ(GR_ENTRY(T, 3 * i, sz)), ctx) == T_TRUE)
            ok = 0;
    if (ok)
    {
        ok = (gr_set(prod, PZ(T), ctx) == GR_SUCCESS);
        for (i = 1; i < num && ok; i++)
            ok = (gr_mul(GR_ENTRY(prod, i, sz), GR_ENTRY(prod, i - 1, sz),
                            PZ(GR_ENTRY(T, 3 * i, sz)), ctx) == GR_SUCCESS);
    }
    if (ok)
        ok = (gr_mul(acc, acc, GR_ENTRY(prod, num - 1, sz), ctx) == GR_SUCCESS);
    if (ok)
        ok = (gr_inv(inv, GR_ENTRY(prod, num - 1, sz), ctx) == GR_SUCCESS);
    for (i = num - 1; i >= 0 && ok; i--)
    {
        gr_ptr Pi = GR_ENTRY(T, 3 * i, sz);
        if (i > 0)
        {
            ok = ok && (gr_mul(t, inv, GR_ENTRY(prod, i - 1, sz), ctx) == GR_SUCCESS);
            ok = ok && (gr_mul(inv, inv, PZ(Pi), ctx) == GR_SUCCESS);
        }
        else
            ok = ok && (gr_set(t, inv, ctx) == GR_SUCCESS);
        ok = ok && (gr_sqr(g, t, ctx) == GR_SUCCESS);
        ok = ok && (gr_mul(PX(Pi), PX(Pi), g, ctx) == GR_SUCCESS);
        ok = ok && (gr_mul(g, g, t, ctx) == GR_SUCCESS);
        ok = ok && (gr_mul(PY(Pi), PY(Pi), g, ctx) == GR_SUCCESS);
        ok = ok && (gr_one(PZ(Pi), ctx) == GR_SUCCESS);
    }
    GR_TMP_CLEAR3(inv, t, g, ctx);
    gr_heap_clear_vec(prod, num, ctx);
    return ok;
}

#define NAF_W 4

/*
    R = k P for an affine point P and k >= 0, multiplying acc by all
    quantities whose invertibility was assumed; returns 0 if acc became
    zero or an operation failed (n composite), else 1. Width-4 signed
    sliding window with the odd multiples in affine coordinates.
*/
int
ecpp_point_mul_gr(gr_ptr R, gr_srcptr P, const fmpz_t k, gr_srcptr a, gr_ptr acc,
                                                                gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem, i, bits = fmpz_bits(k), len;
    const slong ntab = 1 << (NAF_W - 2);
    gr_ptr t1, t2, t3, t4, Q, P2, T, Tn;
    signed char * naf;
    fmpz_t kk;
    int ok = 1;

    if (bits == 0 || gr_is_zero(PZ(P), ctx) == T_TRUE)
        return gr_zero(PZ(R), ctx) == GR_SUCCESS;

    GR_TMP_INIT4(t1, t2, t3, t4, ctx);

    if (bits < 32)
    {
        ok = _gpoint_mul_binary(R, P, k, a, acc, t1, t2, t3, t4, ctx);
        GR_TMP_CLEAR4(t1, t2, t3, t4, ctx);
        return ok && gr_is_zero(acc, ctx) != T_TRUE;
    }

    Q = gr_heap_init_vec(3, ctx);
    P2 = gr_heap_init_vec(3, ctx);
    T = gr_heap_init_vec(3 * ntab, ctx);
    Tn = gr_heap_init_vec(3 * ntab, ctx);

    /* odd multiples: 2P affine, T[i] = (2i + 1) P, Tn[i] = -T[i] */
    ok = _gpoint_double(P2, P, a, acc, t1, t2, t3, t4, ctx);
    if (ok && !_gpoints_to_affine(P2, 1, acc, ctx))
        goto fallback;
    if (ok) ok = _gpoint_set(T, P, ctx);
    for (i = 1; i < ntab && ok; i++)
        ok = _gpoint_add_affine(GR_ENTRY(T, 3 * i, sz), GR_ENTRY(T, 3 * (i - 1), sz),
                                P2, a, acc, t1, t2, t3, t4, ctx);
    if (ok && !_gpoints_to_affine(GR_ENTRY(T, 3, sz), ntab - 1, acc, ctx))
        goto fallback;
    for (i = 0; i < ntab && ok; i++)
    {
        gr_ptr Ti = GR_ENTRY(T, 3 * i, sz), Ni = GR_ENTRY(Tn, 3 * i, sz);
        ok = ok && (gr_set(PX(Ni), PX(Ti), ctx) == GR_SUCCESS);
        ok = ok && (gr_neg(PY(Ni), PY(Ti), ctx) == GR_SUCCESS);
        ok = ok && (gr_one(PZ(Ni), ctx) == GR_SUCCESS);
    }
    if (!ok)
        goto cleanup;

    /* NAF digits, least significant first */
    naf = flint_malloc(bits + 2);
    fmpz_init_set(kk, k);
    len = 0;
    while (!fmpz_is_zero(kk))
    {
        if (fmpz_is_odd(kk))
        {
            slong d = fmpz_fdiv_ui(kk, 1 << NAF_W);
            if (d >= (1 << (NAF_W - 1)))
                d -= 1 << NAF_W;
            naf[len++] = (signed char) d;
            fmpz_sub_si(kk, kk, d);
        }
        else
            naf[len++] = 0;
        fmpz_fdiv_q_2exp(kk, kk, 1);
    }
    fmpz_clear(kk);

    ok = (gr_zero(PZ(Q), ctx) == GR_SUCCESS);
    for (i = len - 1; i >= 0 && ok; i--)
    {
        ok = _gpoint_double(Q, Q, a, acc, t1, t2, t3, t4, ctx);
        if (ok && naf[i] > 0)
            ok = _gpoint_add_affine(Q, Q, GR_ENTRY(T, 3 * ((naf[i] - 1) / 2), sz), a, acc, t1, t2, t3, t4, ctx);
        else if (ok && naf[i] < 0)
            ok = _gpoint_add_affine(Q, Q, GR_ENTRY(Tn, 3 * ((-naf[i] - 1) / 2), sz), a, acc, t1, t2, t3, t4, ctx);
        if (gr_is_zero(acc, ctx) == T_TRUE)
            ok = 0;
    }
    flint_free(naf);
    if (ok)
        ok = _gpoint_set(R, Q, ctx);
    goto cleanup;

fallback:
    /* some small multiple of P is O or not normalisable: plain method */
    ok = (gr_is_zero(acc, ctx) != T_TRUE) && _gpoint_mul_binary(R, P, k, a, acc, t1, t2, t3, t4, ctx);

cleanup:
    gr_heap_clear_vec(Q, 3, ctx);
    gr_heap_clear_vec(P2, 3, ctx);
    gr_heap_clear_vec(T, 3 * ntab, ctx);
    gr_heap_clear_vec(Tn, 3 * ntab, ctx);
    GR_TMP_CLEAR4(t1, t2, t3, t4, ctx);
    return ok && gr_is_zero(acc, ctx) != T_TRUE;
}

/*
    The ring for the modulus n: mpn_mod for 2 to MPN_MOD_MAX_LIMBS limbs
    (much faster than fmpz_mod at these sizes), else fmpz_mod. (A context
    wrapping the given fmpz_mod one without copying would be freed by
    gr_ctx_clear; the copy is cheap next to a scalar multiplication.)
*/
void
ecpp_gr_ctx_init(gr_ctx_t gctx, const fmpz_mod_ctx_t ctx)
{
    const fmpz * n = fmpz_mod_ctx_modulus(ctx);
    if (fmpz_size(n) >= MPN_MOD_MIN_LIMBS && fmpz_size(n) <= MPN_MOD_MAX_LIMBS
            && gr_ctx_init_mpn_mod(gctx, n) == GR_SUCCESS)
        return;
    gr_ctx_init_fmpz_mod(gctx, n);
}

/* the fmpz_mod interface of point.c, through gr */
int
ecpp_point_mul(ecpp_point_t R, const ecpp_point_t P, const fmpz_t k,
                        const fmpz_t a, fmpz_t acc, const fmpz_mod_ctx_t ctx)
{
    gr_ctx_t gctx;
    gr_ptr gP, gR, ga, gacc;
    slong sz;
    int ok;

    ecpp_gr_ctx_init(gctx, ctx);
    sz = gctx->sizeof_elem;
    gP = gr_heap_init_vec(3, gctx);
    gR = gr_heap_init_vec(3, gctx);
    ga = gr_heap_init(gctx);
    gacc = gr_heap_init(gctx);
    ok = (gr_set_fmpz(PX(gP), P->X, gctx) == GR_SUCCESS)
      && (gr_set_fmpz(PY(gP), P->Y, gctx) == GR_SUCCESS)
      && (gr_set_fmpz(PZ(gP), P->Z, gctx) == GR_SUCCESS)
      && (gr_set_fmpz(ga, a, gctx) == GR_SUCCESS)
      && (gr_set_fmpz(gacc, acc, gctx) == GR_SUCCESS);
    if (ok)
        ok = ecpp_point_mul_gr(gR, gP, k, ga, gacc, gctx);
    /* acc and R are read back even on failure (acc may be zero) */
    GR_MUST_SUCCEED(gr_get_fmpz(acc, gacc, gctx));
    if (ok)
    {
        GR_MUST_SUCCEED(gr_get_fmpz(R->X, PX(gR), gctx));
        GR_MUST_SUCCEED(gr_get_fmpz(R->Y, PY(gR), gctx));
        GR_MUST_SUCCEED(gr_get_fmpz(R->Z, PZ(gR), gctx));
    }
    else
        fmpz_zero(R->Z);
    gr_heap_clear_vec(gP, 3, gctx);
    gr_heap_clear_vec(gR, 3, gctx);
    gr_heap_clear(ga, gctx);
    gr_heap_clear(gacc, gctx);
    gr_ctx_clear(gctx);
    return ok;
}
