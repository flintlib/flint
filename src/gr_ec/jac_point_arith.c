/*
    Copyright (C) 2026 Maël Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_ec.h"
#include "gr_ec/impl.h"

#define JX(P) GR_EC_JAC_POINT_X(P, ctx)
#define JY(P) GR_EC_JAC_POINT_Y(P, ctx)
#define JZ(P) GR_EC_JAC_POINT_Z(P, ctx)

#define AX(P) GR_EC_AFF_POINT_X(P, ctx)
#define AY(P) GR_EC_AFF_POINT_Y(P, ctx)

#define T(i) GR_ENTRY(t, i, sz)

int
gr_ec_jac_point_neg(gr_ec_jac_point_t res, const gr_ec_jac_point_t P,
        gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    gr_ptr u, v;
    int status = GR_SUCCESS;

    if (P->is_infinity == T_TRUE)
        return gr_ec_jac_point_zero(res, ctx);
    if (P->is_infinity == T_UNKNOWN)
    {
        res->is_infinity = T_UNKNOWN;
        return GR_UNABLE;
    }

    if (gr_ec_ctx_model(ctx) == GR_EC_SHORT_WEIERSTRASS)
    {
        status |= gr_set(JX(res), JX(P), R);
        status |= gr_neg(JY(res), JY(P), R);
        status |= gr_set(JZ(res), JZ(P), R);
        res->is_infinity = (status == GR_SUCCESS) ? T_FALSE : T_UNKNOWN;
        return status;
    }

    GR_TMP_INIT2(u, v, R);

    /* Y3 = -Y - a1 X Z - a3 Z^3 */
    status |= gr_neg(u, JY(P), R);
    status |= gr_mul(v, JX(P), JZ(P), R);
    status |= gr_submul(u, v, GR_EC_A1(ctx), R);
    status |= gr_sqr(v, JZ(P), R);
    status |= gr_mul(v, v, JZ(P), R);
    status |= gr_submul(u, v, GR_EC_A3(ctx), R);

    status |= gr_set(JX(res), JX(P), R);
    status |= gr_set(JZ(res), JZ(P), R);
    status |= gr_set(JY(res), u, R);
    res->is_infinity = (status == GR_SUCCESS) ? T_FALSE : T_UNKNOWN;

    GR_TMP_CLEAR2(u, v, R);

    return status;
}

/*
    Doubling in the general Weierstrass form, in Jacobian coordinates.
    With x = X/Z^2 and y = Y/Z^3 the affine tangent slope is M/(Z D) where

        D = 2Y + a1 X Z + a3 Z^3,
        M = 3X^2 + 2 a2 X Z^2 + a4 Z^4 - a1 Y Z,

    and taking Z3 = Z D clears all denominators. Reduces to the classical
    short Weierstrass doubling when a1 = a2 = a3 = 0. Assumes P is finite.
*/
int
_gr_ec_jac_point_dbl_long_weierstrass_ws(gr_ec_jac_point_t res,
        const gr_ec_jac_point_t P, gr_ptr t, gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    slong sz = R->sizeof_elem;
    truth_t zero;
    int status = GR_SUCCESS;

    /* T0 = ZZ, T1 = ZZZ, T2 = D = 2Y + a1 X Z + a3 Z^3 */
    status |= gr_sqr(T(0), JZ(P), R);
    status |= gr_mul(T(1), T(0), JZ(P), R);
    status |= gr_add(T(2), JY(P), JY(P), R);
    status |= gr_mul(T(3), JX(P), JZ(P), R);
    status |= gr_addmul(T(2), T(3), GR_EC_A1(ctx), R);
    status |= gr_addmul(T(2), T(1), GR_EC_A3(ctx), R);

    if (status != GR_SUCCESS)
    {
        res->is_infinity = T_UNKNOWN;
        return status;
    }

    zero = gr_is_zero(T(2), R);

    if (zero != T_FALSE)
    {

        if (zero == T_TRUE)
            return gr_ec_jac_point_zero(res, ctx);

        res->is_infinity = T_UNKNOWN;
        return GR_UNABLE;
    }

    /* T3 = M = 3 X^2 + 2 a2 X ZZ + a4 ZZ^2 - a1 Y Z */
    status |= gr_sqr(T(4), JX(P), R);
    status |= gr_add(T(3), T(4), T(4), R);
    status |= gr_add(T(3), T(3), T(4), R);
    status |= gr_mul(T(4), JX(P), T(0), R);
    status |= gr_add(T(4), T(4), T(4), R);
    status |= gr_addmul(T(3), T(4), GR_EC_A2(ctx), R);
    status |= gr_sqr(T(4), T(0), R);
    status |= gr_addmul(T(3), T(4), GR_EC_A4(ctx), R);
    status |= gr_mul(T(4), JY(P), JZ(P), R);
    status |= gr_submul(T(3), T(4), GR_EC_A1(ctx), R);

    /* T4 = v = Z D = Z3, T5 = DD */
    status |= gr_mul(T(4), JZ(P), T(2), R);
    status |= gr_sqr(T(5), T(2), R);

    /* T6 = X3 = M^2 + a1 M v - a2 v^2 - 2 X DD */
    status |= gr_sqr(T(6), T(3), R);
    status |= gr_mul(T(7), T(3), T(4), R);
    status |= gr_addmul(T(6), T(7), GR_EC_A1(ctx), R);
    status |= gr_sqr(T(7), T(4), R);
    status |= gr_submul(T(6), T(7), GR_EC_A2(ctx), R);
    status |= gr_mul(T(7), JX(P), T(5), R);
    status |= gr_sub(T(6), T(6), T(7), R);
    status |= gr_sub(T(6), T(6), T(7), R);

    /* T8 = Y3 = -(M + a1 v) X3 - DD (Y D - M X) - a3 v^3 */
    status |= gr_mul(T(7), GR_EC_A1(ctx), T(4), R);
    status |= gr_add(T(7), T(7), T(3), R);
    status |= gr_mul(T(8), T(7), T(6), R);
    status |= gr_neg(T(8), T(8), R);

    status |= gr_mul(T(7), JY(P), T(2), R);
    status |= gr_submul(T(7), T(3), JX(P), R);
    status |= gr_mul(T(7), T(7), T(5), R);
    status |= gr_sub(T(8), T(8), T(7), R);

    status |= gr_sqr(T(7), T(4), R);
    status |= gr_mul(T(7), T(7), T(4), R);
    status |= gr_submul(T(8), T(7), GR_EC_A3(ctx), R);

    status |= gr_set(JX(res), T(6), R);
    status |= gr_set(JY(res), T(8), R);
    status |= gr_set(JZ(res), T(4), R);
    res->is_infinity = (status == GR_SUCCESS) ? T_FALSE : T_UNKNOWN;


    return status;
}

int
_gr_ec_jac_point_dbl_long_weierstrass(gr_ec_jac_point_t res, const gr_ec_jac_point_t P,
        gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    gr_ptr t;
    int status;

    GR_TMP_INIT_VEC(t, 9, R);
    status = _gr_ec_jac_point_dbl_long_weierstrass_ws(res, P, t, ctx);
    GR_TMP_CLEAR_VEC(t, 9, R);

    return status;
}

/* dbl-2001-b style; assumes P finite and a1 = a2 = a3 = 0 */
int
_gr_ec_jac_point_dbl_short_weierstrass_ws(gr_ec_jac_point_t res,
        const gr_ec_jac_point_t P, gr_ptr t, gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    slong sz = R->sizeof_elem;
    truth_t zero;
    int status = GR_SUCCESS;

    /*
        Everything that reads P is done before the first write to res, so
        that res may alias P and the results can go straight into res
        instead of through a third set of temporaries.
    */

    /* T1 = M = 3 X^2 + a4 Z^4 */
    status |= gr_sqr(T(2), JZ(P), R);
    status |= gr_sqr(T(2), T(2), R);
    status |= gr_sqr(T(0), JX(P), R);
    status |= gr_add(T(1), T(0), T(0), R);
    status |= gr_add(T(1), T(1), T(0), R);
    status |= gr_addmul(T(1), T(2), GR_EC_A4(ctx), R);

    /* T3 = S = 4 X Y^2, T2 = 8 Y^4 */
    status |= gr_sqr(T(2), JY(P), R);
    status |= gr_mul(T(3), JX(P), T(2), R);
    status |= gr_add(T(3), T(3), T(3), R);
    status |= gr_add(T(3), T(3), T(3), R);
    status |= gr_sqr(T(2), T(2), R);
    status |= gr_add(T(2), T(2), T(2), R);
    status |= gr_add(T(2), T(2), T(2), R);
    status |= gr_add(T(2), T(2), T(2), R);

    /* Z3 = 2 Y Z; the last read of P */
    status |= gr_mul(JZ(res), JY(P), JZ(P), R);
    status |= gr_add(JZ(res), JZ(res), JZ(res), R);

    if (status != GR_SUCCESS)
    {
        res->is_infinity = T_UNKNOWN;
        return status;
    }

    zero = gr_is_zero(JZ(res), R);

    if (zero != T_FALSE)
    {
        if (zero == T_TRUE)
            return gr_ec_jac_point_zero(res, ctx);

        res->is_infinity = T_UNKNOWN;
        return GR_UNABLE;
    }

    /* X3 = M^2 - 2 S */
    status |= gr_sqr(JX(res), T(1), R);
    status |= gr_sub(JX(res), JX(res), T(3), R);
    status |= gr_sub(JX(res), JX(res), T(3), R);

    /* Y3 = M (S - X3) - 8 Y^4 */
    status |= gr_sub(T(3), T(3), JX(res), R);
    status |= gr_mul(T(3), T(3), T(1), R);
    status |= gr_sub(JY(res), T(3), T(2), R);

    res->is_infinity = (status == GR_SUCCESS) ? T_FALSE : T_UNKNOWN;

    return status;
}

int
_gr_ec_jac_point_dbl_short_weierstrass(gr_ec_jac_point_t res, const gr_ec_jac_point_t P,
        gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    gr_ptr t;
    int status;

    GR_TMP_INIT_VEC(t, 5, R);
    status = _gr_ec_jac_point_dbl_short_weierstrass_ws(res, P, t, ctx);
    GR_TMP_CLEAR_VEC(t, 5, R);

    return status;
}

int
gr_ec_jac_point_dbl(gr_ec_jac_point_t res, const gr_ec_jac_point_t P,
        gr_ec_ctx_t ctx)
{
    if (P->is_infinity == T_TRUE)
        return gr_ec_jac_point_zero(res, ctx);
    if (P->is_infinity == T_UNKNOWN)
    {
        res->is_infinity = T_UNKNOWN;
        return GR_UNABLE;
    }

    if (gr_ec_ctx_model(ctx) == GR_EC_SHORT_WEIERSTRASS)
        return _gr_ec_jac_point_dbl_short_weierstrass(res, P, ctx);

    return _gr_ec_jac_point_dbl_long_weierstrass(res, P, ctx);
}

/*
    Addition in the general Weierstrass form, in Jacobian coordinates.
    With H = U2 - U1 and r = S2 - S1 as usual, the affine slope is
    r / (H Z1 Z2), and taking Z3 = H Z1 Z2 clears all denominators.
    Reduces to add-1998-cmo when a1 = a2 = a3 = 0. Assumes P, Q finite.
*/
int
_gr_ec_jac_point_add_long_weierstrass(gr_ec_jac_point_t res,
        const gr_ec_jac_point_t P, const gr_ec_jac_point_t Q, gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    slong sz = R->sizeof_elem;
    gr_ptr t;
    truth_t zero;
    int status = GR_SUCCESS;

    GR_TMP_INIT_VEC(t, 11, R);

    /* T0 = Z1^2, T1 = Z2^2, T2 = U1, T3 = U2, T4 = S1, T5 = S2 */
    status |= gr_sqr(T(0), JZ(P), R);
    status |= gr_sqr(T(1), JZ(Q), R);
    status |= gr_mul(T(2), JX(P), T(1), R);
    status |= gr_mul(T(3), JX(Q), T(0), R);
    status |= gr_mul(T(4), JY(P), JZ(Q), R);
    status |= gr_mul(T(4), T(4), T(1), R);
    status |= gr_mul(T(5), JY(Q), JZ(P), R);
    status |= gr_mul(T(5), T(5), T(0), R);

    /* T6 = H = U2 - U1, T5 = r = S2 - S1 */
    status |= gr_sub(T(6), T(3), T(2), R);
    status |= gr_sub(T(5), T(5), T(4), R);

    if (status != GR_SUCCESS)
    {
        GR_TMP_CLEAR_VEC(t, 11, R);
        res->is_infinity = T_UNKNOWN;
        return status;
    }

    zero = gr_is_zero(T(6), R);

    if (zero != T_FALSE)
    {
        truth_t same;

        if (zero == T_UNKNOWN)
        {
            GR_TMP_CLEAR_VEC(t, 11, R);
            res->is_infinity = T_UNKNOWN;
            return GR_UNABLE;
        }

        same = gr_is_zero(T(5), R);
        GR_TMP_CLEAR_VEC(t, 11, R);

        if (same == T_TRUE)
            return _gr_ec_jac_point_dbl_long_weierstrass(res, P, ctx);
        if (same == T_FALSE)
            return gr_ec_jac_point_zero(res, ctx);

        res->is_infinity = T_UNKNOWN;
        return GR_UNABLE;
    }

    /* T7 = v = H Z1 Z2 = Z3, T8 = HH */
    status |= gr_mul(T(7), T(6), JZ(P), R);
    status |= gr_mul(T(7), T(7), JZ(Q), R);
    status |= gr_sqr(T(8), T(6), R);

    /* T9 = X3 = r^2 + a1 r v - a2 v^2 - HH (U1 + U2) */
    status |= gr_sqr(T(9), T(5), R);
    status |= gr_mul(T(10), T(5), T(7), R);
    status |= gr_addmul(T(9), T(10), GR_EC_A1(ctx), R);
    status |= gr_sqr(T(10), T(7), R);
    status |= gr_submul(T(9), T(10), GR_EC_A2(ctx), R);
    status |= gr_add(T(10), T(2), T(3), R);
    status |= gr_mul(T(10), T(10), T(8), R);
    status |= gr_sub(T(9), T(9), T(10), R);

    /* T3 = Y3 = -(r + a1 v) X3 - HH (H S1 - r U1) - a3 v^3 */
    status |= gr_mul(T(10), GR_EC_A1(ctx), T(7), R);
    status |= gr_add(T(10), T(10), T(5), R);
    status |= gr_mul(T(3), T(10), T(9), R);
    status |= gr_neg(T(3), T(3), R);

    status |= gr_mul(T(10), T(6), T(4), R);
    status |= gr_submul(T(10), T(5), T(2), R);
    status |= gr_mul(T(10), T(10), T(8), R);
    status |= gr_sub(T(3), T(3), T(10), R);

    status |= gr_sqr(T(10), T(7), R);
    status |= gr_mul(T(10), T(10), T(7), R);
    status |= gr_submul(T(3), T(10), GR_EC_A3(ctx), R);

    status |= gr_set(JX(res), T(9), R);
    status |= gr_set(JY(res), T(3), R);
    status |= gr_set(JZ(res), T(7), R);
    res->is_infinity = (status == GR_SUCCESS) ? T_FALSE : T_UNKNOWN;

    GR_TMP_CLEAR_VEC(t, 11, R);

    return status;
}

/* add-1998-cmo; assumes P, Q finite and a1 = a2 = a3 = 0 */
int
_gr_ec_jac_point_add_short_weierstrass(gr_ec_jac_point_t res,
        const gr_ec_jac_point_t P, const gr_ec_jac_point_t Q, gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    slong sz = R->sizeof_elem;
    gr_ptr t;
    truth_t zero;
    int status = GR_SUCCESS;

    GR_TMP_INIT_VEC(t, 9, R);

    /* T0 = Z1^2, T1 = Z2^2, T2 = U1, T3 = U2, T4 = S1, T5 = S2 */
    status |= gr_sqr(T(0), JZ(P), R);
    status |= gr_sqr(T(1), JZ(Q), R);
    status |= gr_mul(T(2), JX(P), T(1), R);
    status |= gr_mul(T(3), JX(Q), T(0), R);
    status |= gr_mul(T(4), JY(P), JZ(Q), R);
    status |= gr_mul(T(4), T(4), T(1), R);
    status |= gr_mul(T(5), JY(Q), JZ(P), R);
    status |= gr_mul(T(5), T(5), T(0), R);

    /* T3 = H = U2 - U1, T5 = r = S2 - S1 */
    status |= gr_sub(T(3), T(3), T(2), R);
    status |= gr_sub(T(5), T(5), T(4), R);

    if (status != GR_SUCCESS)
    {
        GR_TMP_CLEAR_VEC(t, 9, R);
        res->is_infinity = T_UNKNOWN;
        return status;
    }

    zero = gr_is_zero(T(3), R);

    if (zero != T_FALSE)
    {
        truth_t same;

        if (zero == T_UNKNOWN)
        {
            GR_TMP_CLEAR_VEC(t, 9, R);
            res->is_infinity = T_UNKNOWN;
            return GR_UNABLE;
        }

        same = gr_is_zero(T(5), R);
        GR_TMP_CLEAR_VEC(t, 9, R);

        if (same == T_TRUE)
            return _gr_ec_jac_point_dbl_short_weierstrass(res, P, ctx);
        if (same == T_FALSE)
            return gr_ec_jac_point_zero(res, ctx);

        res->is_infinity = T_UNKNOWN;
        return GR_UNABLE;
    }

    /* T6 = HH, T7 = HHH, T2 = V = U1 HH */
    status |= gr_sqr(T(6), T(3), R);
    status |= gr_mul(T(7), T(6), T(3), R);
    status |= gr_mul(T(2), T(2), T(6), R);

    /* T8 = X3 = r^2 - HHH - 2 V */
    status |= gr_sqr(T(8), T(5), R);
    status |= gr_sub(T(8), T(8), T(7), R);
    status |= gr_sub(T(8), T(8), T(2), R);
    status |= gr_sub(T(8), T(8), T(2), R);

    /* Y3 = r (V - X3) - S1 HHH */
    status |= gr_sub(T(2), T(2), T(8), R);
    status |= gr_mul(T(2), T(2), T(5), R);
    status |= gr_mul(T(4), T(4), T(7), R);
    status |= gr_sub(T(2), T(2), T(4), R);

    /* Z3 = Z1 Z2 H */
    status |= gr_mul(T(0), JZ(P), JZ(Q), R);
    status |= gr_mul(T(0), T(0), T(3), R);

    status |= gr_set(JX(res), T(8), R);
    status |= gr_set(JY(res), T(2), R);
    status |= gr_set(JZ(res), T(0), R);
    res->is_infinity = (status == GR_SUCCESS) ? T_FALSE : T_UNKNOWN;

    GR_TMP_CLEAR_VEC(t, 9, R);

    return status;
}

int
gr_ec_jac_point_add(gr_ec_jac_point_t res, const gr_ec_jac_point_t P,
        const gr_ec_jac_point_t Q, gr_ec_ctx_t ctx)
{
    if (P->is_infinity == T_TRUE)
        return gr_ec_jac_point_set(res, Q, ctx);
    if (Q->is_infinity == T_TRUE)
        return gr_ec_jac_point_set(res, P, ctx);
    if (P->is_infinity == T_UNKNOWN || Q->is_infinity == T_UNKNOWN)
    {
        res->is_infinity = T_UNKNOWN;
        return GR_UNABLE;
    }

    if (gr_ec_ctx_model(ctx) == GR_EC_SHORT_WEIERSTRASS)
        return _gr_ec_jac_point_add_short_weierstrass(res, P, Q, ctx);

    return _gr_ec_jac_point_add_long_weierstrass(res, P, Q, ctx);
}

int
gr_ec_jac_point_sub(gr_ec_jac_point_t res, const gr_ec_jac_point_t P,
        const gr_ec_jac_point_t Q, gr_ec_ctx_t ctx)
{
    gr_ec_jac_point_t T;
    int status;

    gr_ec_jac_point_init(T, ctx);

    status = gr_ec_jac_point_neg(T, Q, ctx);

    if (status == GR_SUCCESS)
        status = gr_ec_jac_point_add(res, P, T, ctx);

    gr_ec_jac_point_clear(T, ctx);

    return status;
}

/*
    Mixed addition in the general Weierstrass form: the Z2 = 1 case of
    _gr_ec_jac_point_add_long_weierstrass. Assumes P and Q finite.
*/
int
_gr_ec_jac_point_add_aff_point_long_weierstrass_ws(gr_ec_jac_point_t res,
        const gr_ec_jac_point_t P, const gr_ec_aff_point_t Q, gr_ptr t, gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    slong sz = R->sizeof_elem;
    truth_t zero;
    int status = GR_SUCCESS;

    /* T0 = Z1^2, T1 = U2 = X2 Z1^2, T2 = r = Y2 Z1^3 - Y1, T3 = H = U2 - X1 */
    status |= gr_sqr(T(0), JZ(P), R);
    status |= gr_mul(T(1), AX(Q), T(0), R);
    status |= gr_mul(T(2), T(0), JZ(P), R);
    status |= gr_mul(T(2), T(2), AY(Q), R);
    status |= gr_sub(T(2), T(2), JY(P), R);
    status |= gr_sub(T(3), T(1), JX(P), R);

    if (status != GR_SUCCESS)
    {
        res->is_infinity = T_UNKNOWN;
        return status;
    }

    zero = gr_is_zero(T(3), R);

    if (zero != T_FALSE)
    {
        truth_t same;

        if (zero == T_UNKNOWN)
        {
            res->is_infinity = T_UNKNOWN;
            return GR_UNABLE;
        }

        same = gr_is_zero(T(2), R);

        if (same == T_TRUE)
            return _gr_ec_jac_point_dbl_long_weierstrass_ws(res, P, t, ctx);
        if (same == T_FALSE)
            return gr_ec_jac_point_zero(res, ctx);

        res->is_infinity = T_UNKNOWN;
        return GR_UNABLE;
    }

    /* T4 = v = H Z1 = Z3, T5 = HH */
    status |= gr_mul(T(4), T(3), JZ(P), R);
    status |= gr_sqr(T(5), T(3), R);

    /* T6 = X3 = r^2 + a1 r v - a2 v^2 - HH (X1 + U2) */
    status |= gr_sqr(T(6), T(2), R);
    status |= gr_mul(T(7), T(2), T(4), R);
    status |= gr_addmul(T(6), T(7), GR_EC_A1(ctx), R);
    status |= gr_sqr(T(7), T(4), R);
    status |= gr_submul(T(6), T(7), GR_EC_A2(ctx), R);
    status |= gr_add(T(7), JX(P), T(1), R);
    status |= gr_mul(T(7), T(7), T(5), R);
    status |= gr_sub(T(6), T(6), T(7), R);

    /* T8 = Y3 = -(r + a1 v) X3 - HH (H Y1 - r X1) - a3 v^3 */
    status |= gr_mul(T(7), GR_EC_A1(ctx), T(4), R);
    status |= gr_add(T(7), T(7), T(2), R);
    status |= gr_mul(T(8), T(7), T(6), R);
    status |= gr_neg(T(8), T(8), R);

    status |= gr_mul(T(7), T(3), JY(P), R);
    status |= gr_submul(T(7), T(2), JX(P), R);
    status |= gr_mul(T(7), T(7), T(5), R);
    status |= gr_sub(T(8), T(8), T(7), R);

    status |= gr_sqr(T(7), T(4), R);
    status |= gr_mul(T(7), T(7), T(4), R);
    status |= gr_submul(T(8), T(7), GR_EC_A3(ctx), R);

    status |= gr_set(JX(res), T(6), R);
    status |= gr_set(JY(res), T(8), R);
    status |= gr_set(JZ(res), T(4), R);
    res->is_infinity = (status == GR_SUCCESS) ? T_FALSE : T_UNKNOWN;


    return status;
}

int
_gr_ec_jac_point_add_aff_point_long_weierstrass(gr_ec_jac_point_t res, const gr_ec_jac_point_t P,
        const gr_ec_aff_point_t Q,
        gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    gr_ptr t;
    int status;

    GR_TMP_INIT_VEC(t, 9, R);
    status = _gr_ec_jac_point_add_aff_point_long_weierstrass_ws(res, P, Q, t, ctx);
    GR_TMP_CLEAR_VEC(t, 9, R);

    return status;
}

/* madd-2007-bl style; assumes P finite, Q finite and a1 = a2 = a3 = 0 */
int
_gr_ec_jac_point_add_aff_point_short_weierstrass_ws(gr_ec_jac_point_t res,
        const gr_ec_jac_point_t P, const gr_ec_aff_point_t Q, gr_ptr t,
        gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    slong sz = R->sizeof_elem;
    truth_t zero;
    int status = GR_SUCCESS;

    /* T0 = Z1^2, T1 = H = X2 Z1^2 - X1, T2 = r = Y2 Z1^3 - Y1 */
    status |= gr_sqr(T(0), JZ(P), R);
    status |= gr_mul(T(1), AX(Q), T(0), R);
    status |= gr_sub(T(1), T(1), JX(P), R);
    status |= gr_mul(T(2), T(0), JZ(P), R);
    status |= gr_mul(T(2), T(2), AY(Q), R);
    status |= gr_sub(T(2), T(2), JY(P), R);

    if (status != GR_SUCCESS)
    {
        res->is_infinity = T_UNKNOWN;
        return status;
    }

    zero = gr_is_zero(T(1), R);

    if (zero != T_FALSE)
    {
        truth_t same;

        if (zero == T_UNKNOWN)
        {
            res->is_infinity = T_UNKNOWN;
            return GR_UNABLE;
        }

        same = gr_is_zero(T(2), R);

        if (same == T_TRUE)
            return _gr_ec_jac_point_dbl_short_weierstrass_ws(res, P, t, ctx);
        if (same == T_FALSE)
            return gr_ec_jac_point_zero(res, ctx);

        res->is_infinity = T_UNKNOWN;
        return GR_UNABLE;
    }

    /* T4 = HH, T5 = HHH, T6 = HHH Y1, T4 = V = X1 HH */
    status |= gr_sqr(T(4), T(1), R);
    status |= gr_mul(T(5), T(4), T(1), R);
    status |= gr_mul(T(6), T(5), JY(P), R);
    status |= gr_mul(T(4), T(4), JX(P), R);

    /* Z3 = Z1 H; the last read of P */
    status |= gr_mul(JZ(res), JZ(P), T(1), R);

    /* X3 = r^2 - HHH - 2 V */
    status |= gr_sqr(JX(res), T(2), R);
    status |= gr_sub(JX(res), JX(res), T(5), R);
    status |= gr_sub(JX(res), JX(res), T(4), R);
    status |= gr_sub(JX(res), JX(res), T(4), R);

    /* Y3 = r (V - X3) - Y1 HHH */
    status |= gr_sub(T(4), T(4), JX(res), R);
    status |= gr_mul(T(4), T(4), T(2), R);
    status |= gr_sub(JY(res), T(4), T(6), R);

    res->is_infinity = (status == GR_SUCCESS) ? T_FALSE : T_UNKNOWN;

    return status;
}

int
_gr_ec_jac_point_add_aff_point_short_weierstrass(gr_ec_jac_point_t res, const gr_ec_jac_point_t P,
        const gr_ec_aff_point_t Q,
        gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    gr_ptr t;
    int status;

    GR_TMP_INIT_VEC(t, 7, R);
    status = _gr_ec_jac_point_add_aff_point_short_weierstrass_ws(res, P, Q, t, ctx);
    GR_TMP_CLEAR_VEC(t, 7, R);

    return status;
}

int
gr_ec_jac_point_add_aff_point(gr_ec_jac_point_t res, const gr_ec_jac_point_t P,
        const gr_ec_aff_point_t Q, gr_ec_ctx_t ctx)
{
    if (Q->is_infinity == T_TRUE)
        return gr_ec_jac_point_set(res, P, ctx);
    if (P->is_infinity == T_TRUE)
        return gr_ec_jac_point_set_aff_point(res, Q, ctx);
    if (P->is_infinity == T_UNKNOWN || Q->is_infinity == T_UNKNOWN)
    {
        res->is_infinity = T_UNKNOWN;
        return GR_UNABLE;
    }

    if (gr_ec_ctx_model(ctx) == GR_EC_SHORT_WEIERSTRASS)
        return _gr_ec_jac_point_add_aff_point_short_weierstrass(res, P, Q, ctx);

    return _gr_ec_jac_point_add_aff_point_long_weierstrass(res, P, Q, ctx);
}

int
gr_ec_jac_point_sub_aff_point(gr_ec_jac_point_t res, const gr_ec_jac_point_t P,
        const gr_ec_aff_point_t Q, gr_ec_ctx_t ctx)
{
    gr_ec_aff_point_t T;
    int status;

    gr_ec_aff_point_init(T, ctx);

    status = gr_ec_aff_point_neg(T, Q, ctx);

    if (status == GR_SUCCESS)
        status = gr_ec_jac_point_add_aff_point(res, P, T, ctx);

    gr_ec_aff_point_clear(T, ctx);

    return status;
}

int
_gr_ec_jac_point_mul_fmpz_binary(gr_ec_jac_point_t res,
        const gr_ec_jac_point_t P, const fmpz_t n, gr_ec_ctx_t ctx)
{
    gr_ec_jac_point_t A, S;
    fmpz_t k;
    slong i, bits;
    int status = GR_SUCCESS;

    if (fmpz_is_zero(n))
        return gr_ec_jac_point_zero(res, ctx);

    gr_ec_jac_point_init(A, ctx);
    gr_ec_jac_point_init(S, ctx);
    fmpz_init(k);
    fmpz_abs(k, n);

    if (fmpz_sgn(n) < 0)
        status |= gr_ec_jac_point_neg(S, P, ctx);
    else
        status |= gr_ec_jac_point_set(S, P, ctx);

    bits = fmpz_bits(k);

    for (i = bits - 1; i >= 0 && status == GR_SUCCESS; i--)
    {
        status |= gr_ec_jac_point_dbl(A, A, ctx);

        if (fmpz_tstbit(k, i))
            status |= gr_ec_jac_point_add(A, A, S, ctx);
    }

    if (status == GR_SUCCESS)
        status = gr_ec_jac_point_set(res, A, ctx);
    else
        res->is_infinity = T_UNKNOWN;

    fmpz_clear(k);
    gr_ec_jac_point_clear(A, ctx);
    gr_ec_jac_point_clear(S, ctx);

    return status;
}

/*
    Scalar multiplication. The width-w NAF ladder over a batch-normalized
    table of odd multiples is the fast path; it needs one inversion, so it
    falls back to the plain binary ladder over base rings where that is not
    available.
*/
int
gr_ec_jac_point_mul_fmpz(gr_ec_jac_point_t res, const gr_ec_jac_point_t P,
        const fmpz_t n, gr_ec_ctx_t ctx)
{
    int status;

    if (fmpz_is_zero(n) || P->is_infinity == T_TRUE)
        return gr_ec_jac_point_zero(res, ctx);

    if (P->is_infinity == T_UNKNOWN)
    {
        res->is_infinity = T_UNKNOWN;
        return GR_UNABLE;
    }

    status = _gr_ec_jac_point_mul_fmpz_naf(res, P, n, ctx);

    if (status != GR_SUCCESS)
        status = _gr_ec_jac_point_mul_fmpz_binary(res, P, n, ctx);

    return status;
}

int
gr_ec_jac_point_mul_ui(gr_ec_jac_point_t res, const gr_ec_jac_point_t P,
        ulong n, gr_ec_ctx_t ctx)
{
    fmpz_t t;
    int status;

    fmpz_init_set_ui(t, n);
    status = gr_ec_jac_point_mul_fmpz(res, P, t, ctx);
    fmpz_clear(t);

    return status;
}

int
gr_ec_jac_point_mul_si(gr_ec_jac_point_t res, const gr_ec_jac_point_t P,
        slong n, gr_ec_ctx_t ctx)
{
    fmpz_t t;
    int status;

    fmpz_init(t);
    fmpz_set_si(t, n);
    status = gr_ec_jac_point_mul_fmpz(res, P, t, ctx);
    fmpz_clear(t);

    return status;
}
