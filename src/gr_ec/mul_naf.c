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

/*
    Montgomery's trick: normalize a whole vector of Jacobian points with a
    single inversion, at the cost of three multiplications per point. The
    points at infinity in the input contribute nothing to the product and
    come out as points at infinity.
*/
int
_gr_ec_jac_point_vec_get_aff_point_vec(gr_ec_aff_point_struct * res,
        const gr_ec_jac_point_struct * P, slong len, gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    slong sz = R->sizeof_elem;
    gr_ptr c, inv, zi, zz;
    slong i;
    int status = GR_SUCCESS;

    if (len <= 0)
        return GR_SUCCESS;

    for (i = 0; i < len; i++)
        if (P[i].is_infinity == T_UNKNOWN)
            return GR_UNABLE;

    c = gr_heap_init_vec(len, R);
    GR_TMP_INIT3(inv, zi, zz, R);

    /* c[i] = product of the Z of the finite points up to index i */
    status |= gr_one(inv, R);

    for (i = 0; i < len; i++)
    {
        if (P[i].is_infinity == T_FALSE)
            status |= gr_mul(inv, inv, GR_EC_JAC_POINT_Z(&P[i], ctx), R);

        status |= gr_set(GR_ENTRY(c, i, sz), inv, R);
    }

    if (status == GR_SUCCESS)
        status = gr_inv(inv, GR_ENTRY(c, len - 1, sz), R);

    for (i = len - 1; i >= 0 && status == GR_SUCCESS; i--)
    {
        if (P[i].is_infinity != T_FALSE)
        {
            res[i].is_infinity = T_TRUE;
            continue;
        }

        if (i > 0)
        {
            status |= gr_mul(zi, inv, GR_ENTRY(c, i - 1, sz), R);
            status |= gr_mul(inv, inv, GR_EC_JAC_POINT_Z(&P[i], ctx), R);
        }
        else
            status |= gr_set(zi, inv, R);

        status |= gr_sqr(zz, zi, R);
        status |= gr_mul(GR_EC_AFF_POINT_X(&res[i], ctx),
                    GR_EC_JAC_POINT_X(&P[i], ctx), zz, R);
        status |= gr_mul(zz, zz, zi, R);
        status |= gr_mul(GR_EC_AFF_POINT_Y(&res[i], ctx),
                    GR_EC_JAC_POINT_Y(&P[i], ctx), zz, R);
        res[i].is_infinity = T_FALSE;
    }

    if (status != GR_SUCCESS)
        for (i = 0; i < len; i++)
            res[i].is_infinity = T_UNKNOWN;

    GR_TMP_CLEAR3(inv, zi, zz, R);
    gr_heap_clear_vec(c, len, R);

    return status;
}

/*
    Width-w non-adjacent form of a positive integer: digits[i] is the
    coefficient of 2^i and is either zero or an odd integer of absolute
    value less than 2^(w-1). Returns the number of digits written, which
    is at most fmpz_bits(n) + 1.
*/
static slong
_gr_ec_wnaf(slong * digits, const fmpz_t n, slong w)
{
    fmpz_t k;
    slong len = 0;
    ulong mod = UWORD(1) << w;
    ulong half = UWORD(1) << (w - 1);

    fmpz_init_set(k, n);

    while (!fmpz_is_zero(k))
    {
        if (fmpz_is_odd(k))
        {
            slong d = (slong) fmpz_fdiv_ui(k, mod);

            if ((ulong) d >= half)
                d -= (slong) mod;

            digits[len] = d;
            fmpz_sub_si(k, k, d);
        }
        else
            digits[len] = 0;

        len++;
        fmpz_fdiv_q_2exp(k, k, 1);
    }

    fmpz_clear(k);

    return len;
}

/* window width; the table holds 2^(w-2) affine points */
static slong
_gr_ec_naf_window(slong bits)
{
    if (bits < 128)
        return 3;
    if (bits < 192)
        return 4;
    if (bits < 1024)
        return 5;
    return 6;
}

/*
    Scalar multiplication by a width-w NAF ladder over a table of odd
    multiples P, 3P, 5P, ..., which is normalized to affine coordinates
    with a single inversion so that every addition in the main loop is a
    mixed addition. Returns GR_UNABLE if the table cannot be normalized,
    which is the case over a base ring that is not a field; the caller
    should fall back to _gr_ec_jac_point_mul_fmpz_binary there.
*/
int
_gr_ec_jac_point_mul_fmpz_naf(gr_ec_jac_point_t res, const gr_ec_jac_point_t P,
        const fmpz_t n, gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    slong sz = R->sizeof_elem;
    gr_ec_jac_point_struct * T;
    gr_ec_aff_point_struct * A;
    gr_ec_jac_point_struct * D;
    gr_ec_jac_point_struct * acc;
    gr_ptr jac_coords, aff_coords, t;
    slong * digits;
    fmpz_t k;
    slong w, ntab, i, len, bits;
    int shrt, status = GR_SUCCESS;

    if (fmpz_is_zero(n) || P->is_infinity == T_TRUE)
        return gr_ec_jac_point_zero(res, ctx);

    if (P->is_infinity == T_UNKNOWN)
    {
        res->is_infinity = T_UNKNOWN;
        return GR_UNABLE;
    }

    fmpz_init(k);
    fmpz_abs(k, n);
    bits = fmpz_bits(k);

    w = _gr_ec_naf_window(bits);
    ntab = WORD(1) << (w - 2);
    shrt = (gr_ec_ctx_model(ctx) == GR_EC_SHORT_WEIERSTRASS);

    /* the table, the accumulator and the workspace in two allocations */
    jac_coords = gr_heap_init_vec(3 * (ntab + 2), R);
    aff_coords = gr_heap_init_vec(2 * (2 * ntab), R);
    T = flint_malloc((ntab + 2) * sizeof(gr_ec_jac_point_struct));
    A = flint_malloc(2 * ntab * sizeof(gr_ec_aff_point_struct));
    digits = flint_malloc((bits + 2) * sizeof(slong));

    for (i = 0; i < ntab + 2; i++)
    {
        T[i].coords = GR_ENTRY(jac_coords, 3 * i, sz);
        T[i].is_infinity = T_TRUE;
    }

    for (i = 0; i < 2 * ntab; i++)
    {
        A[i].coords = GR_ENTRY(aff_coords, 2 * i, sz);
        A[i].is_infinity = T_TRUE;
    }

    D = &T[ntab];
    acc = &T[ntab + 1];

    /* T[i] = (2i + 1) P */
    status |= gr_ec_jac_point_set(&T[0], P, ctx);

    if (ntab > 1)
    {
        status |= gr_ec_jac_point_dbl(D, P, ctx);

        for (i = 1; i < ntab && status == GR_SUCCESS; i++)
            status |= gr_ec_jac_point_add(&T[i], &T[i - 1], D, ctx);
    }

    /* one inversion normalizes the whole table */
    if (status == GR_SUCCESS)
        status = _gr_ec_jac_point_vec_get_aff_point_vec(A, T, ntab, ctx);

    /* and the negated table, so the ladder never negates */
    for (i = 0; i < ntab && status == GR_SUCCESS; i++)
        status |= (fmpz_sgn(n) < 0) ? gr_ec_aff_point_neg(&A[i], &A[i], ctx)
                                    : GR_SUCCESS;

    for (i = 0; i < ntab && status == GR_SUCCESS; i++)
        status |= gr_ec_aff_point_neg(&A[ntab + i], &A[i], ctx);

    if (status == GR_SUCCESS)
    {
        len = _gr_ec_wnaf(digits, k, w);

        GR_TMP_INIT_VEC(t, GR_EC_JAC_SCRATCH, R);

        for (i = len - 1; i >= 0 && status == GR_SUCCESS; i--)
        {
            const gr_ec_aff_point_struct * S;
            slong d = digits[i];

            if (acc->is_infinity == T_FALSE)
                status |= shrt ? _gr_ec_jac_point_dbl_short_weierstrass_ws(acc, acc, t, ctx)
                               : _gr_ec_jac_point_dbl_long_weierstrass_ws(acc, acc, t, ctx);
            else if (acc->is_infinity == T_UNKNOWN)
                status |= GR_UNABLE;

            if (d == 0 || status != GR_SUCCESS)
                continue;

            S = (d > 0) ? &A[(d - 1) / 2] : &A[ntab + (-d - 1) / 2];

            if (S->is_infinity == T_TRUE)
                continue;

            if (acc->is_infinity == T_TRUE)
                status |= gr_ec_jac_point_set_aff_point(acc, S, ctx);
            else
                status |= shrt ? _gr_ec_jac_point_add_aff_point_short_weierstrass_ws(acc, acc, S, t, ctx)
                               : _gr_ec_jac_point_add_aff_point_long_weierstrass_ws(acc, acc, S, t, ctx);
        }

        GR_TMP_CLEAR_VEC(t, GR_EC_JAC_SCRATCH, R);
    }

    if (status == GR_SUCCESS)
        status = gr_ec_jac_point_set(res, acc, ctx);
    else
        res->is_infinity = T_UNKNOWN;

    gr_heap_clear_vec(jac_coords, 3 * (ntab + 2), R);
    gr_heap_clear_vec(aff_coords, 2 * (2 * ntab), R);
    flint_free(T);
    flint_free(A);
    flint_free(digits);
    fmpz_clear(k);

    return status;
}
