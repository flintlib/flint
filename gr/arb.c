/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "arb_poly.h"
#include "arb_mat.h"
#include "arb_fmpz_poly.h"
#include "arb_hypgeom.h"
#include "fmpzi.h"
#include "qqbar.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_poly.h"

typedef struct
{
    slong prec;
}
gr_arb_ctx;

#define ARB_CTX_PREC(ring_ctx) (((gr_arb_ctx *)((ring_ctx)))->prec)

int _gr_arb_ctx_set_real_prec(gr_ctx_t ctx, slong prec)
{
    prec = FLINT_MAX(prec, 2);
    prec = FLINT_MIN(prec, WORD_MAX / 8);

    ARB_CTX_PREC(ctx) = prec;
    return GR_SUCCESS;
}

int _gr_arb_ctx_get_real_prec(slong * res, gr_ctx_t ctx)
{
    *res = ARB_CTX_PREC(ctx);
    return GR_SUCCESS;
}

int
_gr_arb_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    gr_stream_write(out, "Real numbers (arb, prec = ");
    gr_stream_write_si(out, ARB_CTX_PREC(ctx));
    gr_stream_write(out, ")");
    return GR_SUCCESS;
}

void
_gr_arb_init(arb_t x, const gr_ctx_t ctx)
{
    arb_init(x);
}

void
_gr_arb_clear(arb_t x, const gr_ctx_t ctx)
{
    arb_clear(x);
}

void
_gr_arb_swap(arb_t x, arb_t y, const gr_ctx_t ctx)
{
    arb_t t;
    *t = *x;
    *x = *y;
    *y = *t;
}

/* todo: limits */
int
_gr_arb_randtest(arb_t res, flint_rand_t state, const gr_ctx_t ctx)
{
    arb_randtest(res, state, ARB_CTX_PREC(ctx), 10);
    return GR_SUCCESS;
}

/* todo */
int
_gr_arb_write(gr_stream_t out, const arb_t x, const gr_ctx_t ctx)
{
    gr_stream_write_free(out, arb_get_str(x, ARB_CTX_PREC(ctx) * 0.30102999566398 + 1, 0));
    return GR_SUCCESS;
}

int
_gr_arb_zero(arb_t x, const gr_ctx_t ctx)
{
    arb_zero(x);
    return GR_SUCCESS;
}

int
_gr_arb_one(arb_t x, const gr_ctx_t ctx)
{
    arb_one(x);
    return GR_SUCCESS;
}

int
_gr_arb_set_si(arb_t res, slong v, const gr_ctx_t ctx)
{
    arb_set_si(res, v);
    arb_set_round(res, res, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_set_ui(arb_t res, ulong v, const gr_ctx_t ctx)
{
    arb_set_ui(res, v);
    arb_set_round(res, res, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_set_fmpz(arb_t res, const fmpz_t v, const gr_ctx_t ctx)
{
    arb_set_round_fmpz(res, v, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_set_fmpq(arb_t res, const fmpq_t v, const gr_ctx_t ctx)
{
    arb_set_fmpq(res, v, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_set_str(arb_t res, const char * x, const gr_ctx_t ctx)
{
    if (arb_set_str(res, x, ARB_CTX_PREC(ctx)))
        return GR_DOMAIN;

    return GR_SUCCESS;
}

int
_gr_arb_set_d(arb_t res, double x, const gr_ctx_t ctx)
{
    arb_set_d(res, x);
    arb_set_round(res, res, ARB_CTX_PREC(ctx));

    if (!arb_is_finite(res))
        return GR_DOMAIN;

    return GR_SUCCESS;
}

int 
_gr_ca_get_arb_with_prec(arb_t res, gr_srcptr x, gr_ctx_t x_ctx, slong prec);

int
_gr_arb_set_other(arb_t res, gr_srcptr x, gr_ctx_t x_ctx, const gr_ctx_t ctx)
{
    switch (x_ctx->which_ring)
    {
        case GR_CTX_FMPZ:
            return _gr_arb_set_fmpz(res, x, ctx);

        case GR_CTX_FMPQ:
            return _gr_arb_set_fmpq(res, x, ctx);

        case GR_CTX_FMPZI:
            if (!fmpz_is_zero(fmpzi_imagref((const fmpzi_struct *) x)))
                return GR_DOMAIN;
            arb_set_round_fmpz(res, fmpzi_realref((const fmpzi_struct *) x), ARB_CTX_PREC(ctx));
            return GR_SUCCESS;

        case GR_CTX_REAL_ALGEBRAIC_QQBAR:
            qqbar_get_arb(res, x, ARB_CTX_PREC(ctx));
            return GR_SUCCESS;

        case GR_CTX_COMPLEX_ALGEBRAIC_QQBAR:
            if (qqbar_is_real(x))
            {
                qqbar_get_arb(res, x, ARB_CTX_PREC(ctx));
                return GR_SUCCESS;
            }
            else
            {
                return GR_DOMAIN;
            }

        case GR_CTX_RR_CA:
        case GR_CTX_REAL_ALGEBRAIC_CA:
        case GR_CTX_CC_CA:
        case GR_CTX_COMPLEX_ALGEBRAIC_CA:
            return _gr_ca_get_arb_with_prec(res, x, x_ctx, ARB_CTX_PREC(ctx));

        case GR_CTX_REAL_FLOAT_ARF:
            if (arf_is_finite(x))
            {
                arb_set_arf(res, x);
                arb_set_round(res, res, ARB_CTX_PREC(ctx));
                return GR_SUCCESS;
            }
            else
            {
                return GR_DOMAIN;
            }

        case GR_CTX_RR_ARB:
            arb_set_round(res, x, ARB_CTX_PREC(ctx));
            return GR_SUCCESS;

        case GR_CTX_CC_ACB:
            if (arb_is_zero(acb_imagref((acb_srcptr) x)))
            {
                arb_set_round(res, x, ARB_CTX_PREC(ctx));
                return GR_SUCCESS;
            }
            else if (arb_contains_zero(acb_imagref((acb_srcptr) x)))
            {
                return GR_UNABLE;
            }
            else
            {
                return GR_DOMAIN;
            }
    }

    return GR_UNABLE;
}

/* xxx: assumes that ctx are not read */
int _gr_arf_get_fmpz(fmpz_t res, const arf_t x, const gr_ctx_t ctx);
int _gr_arf_get_si(slong * res, const arf_t x, const gr_ctx_t ctx);
int _gr_arf_get_ui(ulong * res, const arf_t x, const gr_ctx_t ctx);

int
_gr_arb_get_fmpz(fmpz_t res, const arb_t x, const gr_ctx_t ctx)
{
    if (!arb_is_int(x))
    {
        if (arb_contains_int(x))
            return GR_UNABLE;
        else
            return GR_DOMAIN;
    }

    return _gr_arf_get_fmpz(res, arb_midref(x), NULL);
}

int
_gr_arb_get_si(slong * res, const arb_t x, const gr_ctx_t ctx)
{
    if (!arb_is_int(x))
    {
        if (arb_contains_int(x))
            return GR_UNABLE;
        else
            return GR_DOMAIN;
    }

    return _gr_arf_get_si(res, arb_midref(x), NULL);
}

int
_gr_arb_get_ui(ulong * res, const arb_t x, const gr_ctx_t ctx)
{
    if (!arb_is_int(x))
    {
        if (arb_contains_int(x))
            return GR_UNABLE;
        else
            return GR_DOMAIN;
    }

    return _gr_arf_get_ui(res, arb_midref(x), NULL);
}

int
_gr_arb_get_d(double * res, const arb_t x, const gr_ctx_t ctx)
{
    *res = arf_get_d(arb_midref(x), ARF_RND_NEAR);
    return GR_SUCCESS;
}

truth_t
_gr_arb_is_zero(const arb_t x, const gr_ctx_t ctx)
{
    if (arb_is_zero(x))
        return T_TRUE;

    if (arb_contains_zero(x))
        return T_UNKNOWN;

    return T_FALSE;
}

truth_t
_gr_arb_is_one(const arb_t x, const gr_ctx_t ctx)
{
    if (arb_is_one(x))
        return T_TRUE;

    if (arb_contains_si(x, 1))
        return T_UNKNOWN;

    return T_FALSE;
}

truth_t
_gr_arb_is_neg_one(const arb_t x, const gr_ctx_t ctx)
{
    if (arb_equal_si(x, -1))
        return T_TRUE;

    if (arb_contains_si(x, -1))
        return T_UNKNOWN;

    return T_FALSE;
}

truth_t
_gr_arb_equal(const arb_t x, const arb_t y, const gr_ctx_t ctx)
{
    if (arb_is_exact(x) && arb_equal(x, y))
        return T_TRUE;

    if (arb_overlaps(x, y))
        return T_UNKNOWN;

    return T_FALSE;
}

int
_gr_arb_set(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    arb_set(res, x);
    return GR_SUCCESS;
}

int
_gr_arb_neg(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    arb_neg(res, x);
    return GR_SUCCESS;
}

int
_gr_arb_add(arb_t res, const arb_t x, const arb_t y, const gr_ctx_t ctx)
{
    arb_add(res, x, y, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_add_si(arb_t res, const arb_t x, slong y, const gr_ctx_t ctx)
{
    arb_add_si(res, x, y, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_add_ui(arb_t res, const arb_t x, ulong y, const gr_ctx_t ctx)
{
    arb_add_ui(res, x, y, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_add_fmpz(arb_t res, const arb_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    arb_add_fmpz(res, x, y, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_sub(arb_t res, const arb_t x, const arb_t y, const gr_ctx_t ctx)
{
    arb_sub(res, x, y, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_sub_si(arb_t res, const arb_t x, slong y, const gr_ctx_t ctx)
{
    arb_sub_si(res, x, y, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_sub_ui(arb_t res, const arb_t x, ulong y, const gr_ctx_t ctx)
{
    arb_sub_ui(res, x, y, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_sub_fmpz(arb_t res, const arb_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    arb_sub_fmpz(res, x, y, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_mul(arb_t res, const arb_t x, const arb_t y, const gr_ctx_t ctx)
{
    arb_mul(res, x, y, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_mul_si(arb_t res, const arb_t x, slong y, const gr_ctx_t ctx)
{
    arb_mul_si(res, x, y, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_mul_ui(arb_t res, const arb_t x, ulong y, const gr_ctx_t ctx)
{
    arb_mul_ui(res, x, y, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_mul_fmpz(arb_t res, const arb_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    arb_mul_fmpz(res, x, y, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_addmul(arb_t res, const arb_t x, const arb_t y, const gr_ctx_t ctx)
{
    arb_addmul(res, x, y, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_submul(arb_t res, const arb_t x, const arb_t y, const gr_ctx_t ctx)
{
    arb_submul(res, x, y, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_mul_two(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    arb_mul_2exp_si(res, x, 1);
    return GR_SUCCESS;
}

int
_gr_arb_sqr(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    arb_sqr(res, x, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_mul_2exp_si(arb_t res, const arb_t x, slong y, const gr_ctx_t ctx)
{
    arb_mul_2exp_si(res, x, y);
    return GR_SUCCESS;
}

int
_gr_arb_mul_2exp_fmpz(arb_t res, const arb_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    arb_mul_2exp_fmpz(res, x, y);
    return GR_SUCCESS;
}

int
_gr_arb_set_fmpz_2exp_fmpz(arb_t res, const fmpz_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    arb_set_fmpz_2exp(res, x, y);
    return GR_SUCCESS;
}

int
_gr_arb_get_fmpz_2exp_fmpz(fmpz_t res1, fmpz_t res2, const arb_t x, const gr_ctx_t ctx)
{
    if (!arb_is_exact(x))
        return GR_UNABLE;

    if (!arb_is_finite(x))
        return GR_DOMAIN;

    arf_get_fmpz_2exp(res1, res2, arb_midref(x));
    return GR_SUCCESS;
}


int
_gr_arb_inv(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    if (arb_is_zero(x))
    {
        return GR_DOMAIN;
    }
    else
    {
        arb_inv(res, x, ARB_CTX_PREC(ctx));
        if (arb_is_finite(res))
            return GR_SUCCESS;
        else
            return GR_UNABLE;
    }
}

int
_gr_arb_div(arb_t res, const arb_t x, const arb_t y, const gr_ctx_t ctx)
{
    if (arb_is_zero(y))
    {
        return GR_DOMAIN;
    }
    else
    {
        arb_div(res, x, y, ARB_CTX_PREC(ctx));

        if (arb_is_finite(res))
            return GR_SUCCESS;
        else
            return GR_UNABLE;
    }
}

int
_gr_arb_div_si(arb_t res, const arb_t x, slong y, const gr_ctx_t ctx)
{
    if (y == 0)
    {
        return GR_DOMAIN;
    }
    else
    {
        arb_div_si(res, x, y, ARB_CTX_PREC(ctx));
        return GR_SUCCESS;
    }
}

int
_gr_arb_div_ui(arb_t res, const arb_t x, ulong y, const gr_ctx_t ctx)
{
    if (y == 0)
    {
        return GR_DOMAIN;
    }
    else
    {
        arb_div_ui(res, x, y, ARB_CTX_PREC(ctx));
        return GR_SUCCESS;
    }
}

int
_gr_arb_div_fmpz(arb_t res, const arb_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    if (fmpz_is_zero(y))
    {
        return GR_DOMAIN;
    }
    else
    {
        arb_div_fmpz(res, x, y, ARB_CTX_PREC(ctx));
        return GR_SUCCESS;
    }
}

truth_t
_gr_arb_is_invertible(const arb_t x, const gr_ctx_t ctx)
{
    if (arb_is_zero(x))
        return T_FALSE;

    if (arb_contains_zero(x))
        return T_UNKNOWN;

    return T_TRUE;
}

int
_gr_arb_pow_ui(arb_t res, const arb_t x, ulong exp, const gr_ctx_t ctx)
{
    arb_pow_ui(res, x, exp, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_pow_si(arb_t res, const arb_t x, slong exp, const gr_ctx_t ctx)
{
    if (exp < 0 && arb_is_zero(x))
    {
        return GR_DOMAIN;
    }
    else if (exp < 0 && arb_contains_zero(x))
    {
        return GR_UNABLE;
    }
    else
    {
        fmpz_t t;
        fmpz_init_set_si(t, exp);
        arb_pow_fmpz(res, x, t, ARB_CTX_PREC(ctx));
        fmpz_clear(t);
        return GR_SUCCESS;
    }
}

int
_gr_arb_pow_fmpz(arb_t res, const arb_t x, const fmpz_t exp, const gr_ctx_t ctx)
{
    if (fmpz_sgn(exp) < 0 && arb_is_zero(x))
    {
        return GR_DOMAIN;
    }
    else if (fmpz_sgn(exp) < 0 && arb_contains_zero(x))
    {
        return GR_UNABLE;
    }
    else
    {
        arb_pow_fmpz(res, x, exp, ARB_CTX_PREC(ctx));
        return GR_SUCCESS;
    }
}

int
_gr_arb_pow_fmpq(arb_t res, const arb_t x, const fmpq_t exp, const gr_ctx_t ctx)
{
    if (fmpq_sgn(exp) < 0 && arb_is_zero(x))
    {
        return GR_DOMAIN;
    }
    else if (fmpq_sgn(exp) < 0 && arb_contains_zero(x))
    {
        return GR_UNABLE;
    }
    else
    {
        if (fmpz_is_one(fmpq_denref(exp)) || arb_is_nonnegative(x))
        {
            arb_pow_fmpq(res, x, exp, ARB_CTX_PREC(ctx));
            return GR_SUCCESS;
        }
        else if (arb_is_negative(x))
        {
            return GR_DOMAIN;
        }
        else
        {
            return GR_UNABLE;
        }
    }
}

int
_gr_arb_pow(arb_t res, const arb_t x, const arb_t exp, const gr_ctx_t ctx)
{
    if (arb_is_int(exp))
    {
        if (arf_sgn(arb_midref(exp)) < 0)
        {
            if (arb_is_zero(x))
                return GR_DOMAIN;

            if (arb_contains_zero(x))
                return GR_UNABLE;
        }

        arb_pow(res, x, exp, ARB_CTX_PREC(ctx));
        return GR_SUCCESS;
    }
    else if (arb_is_positive(x) || (arb_is_nonnegative(x) && arb_is_nonnegative(exp)))
    {
        arb_pow(res, x, exp, ARB_CTX_PREC(ctx));
        return GR_SUCCESS;
    }
    else if (arb_is_zero(x) && arb_is_negative(exp))
    {
        return GR_DOMAIN;
    }
    else if (arb_is_negative(x) && !arb_contains_int(exp))
    {
        return GR_DOMAIN;
    }
    else
    {
        return GR_UNABLE;
    }
}

truth_t
_gr_arb_is_square(const arb_t x, const gr_ctx_t ctx)
{
    return T_TRUE;
}

int
_gr_arb_sqrt(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    if (arb_is_nonnegative(x))
    {
        arb_sqrt(res, x, ARB_CTX_PREC(ctx));
        return GR_SUCCESS;
    }
    else if (arb_is_negative(x))
    {
        return GR_DOMAIN;
    }
    else
    {
        return GR_UNABLE;
    }
}

int
_gr_arb_rsqrt(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    if (arb_is_positive(x))
    {
        arb_rsqrt(res, x, ARB_CTX_PREC(ctx));
        return GR_SUCCESS;
    }
    else if (arb_is_nonpositive(x))
    {
        return GR_DOMAIN;
    }
    else
    {
        return GR_UNABLE;
    }
}

int
_gr_arb_floor(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    arb_floor(res, x, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_ceil(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    arb_ceil(res, x, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_arb_trunc(arb_t res, const arb_t x, slong prec)
{
    if (arb_contains_zero(x))
    {
        arb_t a;
        arb_init(a);

        mag_one(arb_radref(a));

        if (arb_contains(a, x))
        {
            arb_zero(res);
        }
        else
        {
            arb_t b;
            arb_init(b);
            arb_floor(a, x, prec);
            arb_ceil(b, x, prec);
            arb_union(res, a, b, prec);
            arb_clear(b);
        }

        arb_clear(a);
    }
    else if (arf_sgn(arb_midref(x)) > 0)
    {
        arb_floor(res, x, prec);
    }
    else
    {
        arb_ceil(res, x, prec);
    }

    return GR_SUCCESS;
}

int
_arb_nint(arb_t res, const arb_t x, slong prec)
{
    if (arb_is_int(x))
    {
        arb_set(res, x);
    }
    else
    {
        arb_t t, u;
        arb_init(t);
        arb_init(u);

        arb_set_d(t, 0.5);
        arb_add(t, x, t, prec);

        arb_mul_2exp_si(u, x, 1);
        arb_sub_ui(u, u, 1, prec);
        arb_mul_2exp_si(u, u, -2);

        arb_floor(res, t, prec);

        /* nint(x) = floor(x+0.5) - isint((2*x-1)/4) */

        if (arb_is_int(u))
        {
            arb_sub_ui(res, res, 1, prec);
        }
        else if (arb_contains_int(u))
        {
            arf_one(arb_midref(u));
            mag_one(arb_radref(u));
            arb_mul_2exp_si(u, u, -1);
            arb_sub_ui(res, res, 1, prec);
        }
    }

    return GR_SUCCESS;
}

int
_gr_arb_trunc(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    return _arb_trunc(res, x, ARB_CTX_PREC(ctx));
}

int
_gr_arb_nint(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    return _arb_nint(res, x, ARB_CTX_PREC(ctx));
}

int
_gr_arb_abs(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    arb_abs(res, x);
    return GR_SUCCESS;
}

int
_gr_arb_conj(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    arb_set(res, x);
    return GR_SUCCESS;
}

int
_gr_arb_im(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    arb_zero(res);
    return GR_SUCCESS;
}

int
_gr_arb_sgn(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    arb_sgn(res, x);
    return GR_SUCCESS;
}

int
_gr_arb_cmp(int * res, const arb_t x, const arb_t y, const gr_ctx_t ctx)
{
    if ((arb_is_exact(x) && arb_is_exact(y)) || !arb_overlaps(x, y))
    {
        *res = arf_cmp(arb_midref(x), arb_midref(y));
        return GR_SUCCESS;
    }
    else
    {
        *res = 0;
        return GR_UNABLE;
    }
}

int
_gr_arb_cmpabs(int * res, const arb_t x, const arb_t y, const gr_ctx_t ctx)
{
    arb_t t, u;

    *t = *x;
    *u = *y;

    if (arf_sgn(arb_midref(t)) < 0)
        ARF_NEG(arb_midref(t));

    if (arf_sgn(arb_midref(u)) < 0)
        ARF_NEG(arb_midref(u));

    return _gr_arb_cmp(res, t, u, ctx);
}

int
_gr_arb_pi(arb_t res, const gr_ctx_t ctx)
{
    arb_const_pi(res, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_euler(arb_t res, const gr_ctx_t ctx)
{
    arb_const_euler(res, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_catalan(arb_t res, const gr_ctx_t ctx)
{
    arb_const_catalan(res, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_khinchin(arb_t res, const gr_ctx_t ctx)
{
    arb_const_khinchin(res, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_glaisher(arb_t res, const gr_ctx_t ctx)
{
    arb_const_glaisher(res, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_exp(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    arb_exp(res, x, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_log(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    if (arb_is_positive(x))
    {
        arb_log(res, x, ARB_CTX_PREC(ctx));
        return GR_SUCCESS;
    }

    if (arb_is_nonpositive(x))
        return GR_DOMAIN;

    return GR_UNABLE;
}

int _gr_arb_sin(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    arb_sin(res, x, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_cos(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    arb_cos(res, x, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_tan(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    arb_tan(res, x, ARB_CTX_PREC(ctx));
    return arb_is_finite(res) ? GR_SUCCESS : GR_DOMAIN;
}

int
_gr_arb_atan(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    arb_atan(res, x, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

#include "bernoulli.h"

/* todo: caching, e.g. if (x <= 1000) bernoulli_cache_compute(x + 1); */

int
_gr_arb_bernoulli_ui(arb_t res, ulong x, const gr_ctx_t ctx)
{
    arb_bernoulli_ui(res, x, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_bernoulli_fmpz(arb_t res, const fmpz_t x, const gr_ctx_t ctx)
{
    arb_bernoulli_fmpz(res, x, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_erf(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    arb_hypgeom_erf(res, x, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_erfc(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    arb_hypgeom_erfc(res, x, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_erfi(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    arb_hypgeom_erfi(res, x, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_gamma(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    if (arb_is_positive(x))
    {
        arb_gamma(res, x, ARB_CTX_PREC(ctx));
        return GR_SUCCESS;
    }
    else if (arb_is_nonpositive(x) && arb_is_int(x))
    {
        return GR_DOMAIN;
    }
    else
    {
        arb_gamma(res, x, ARB_CTX_PREC(ctx));
        return arb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
    }
}

int
_gr_arb_rgamma(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    arb_rgamma(res, x, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_lgamma(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    if (arb_is_positive(x))
    {
        arb_lgamma(res, x, ARB_CTX_PREC(ctx));
        return GR_SUCCESS;
    }

    if (arb_is_nonpositive(x))
        return GR_DOMAIN;

    return GR_UNABLE;
}

int
_gr_arb_digamma(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    if (arb_is_positive(x))
    {
        arb_digamma(res, x, ARB_CTX_PREC(ctx));
        return GR_SUCCESS;
    }
    else if (arb_is_nonpositive(x) && arb_is_int(x))
    {
        return GR_DOMAIN;
    }
    else
    {
        arb_digamma(res, x, ARB_CTX_PREC(ctx));
        return arb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
    }
}

int
_gr_arb_zeta(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    if (arb_contains_si(x, 1))
    {
        if (arb_is_one(x))
            return GR_DOMAIN;
        else
            return GR_UNABLE;
    }
    else
    {
        arb_zeta(res, x, ARB_CTX_PREC(ctx));
        return GR_SUCCESS;
    }
}

int
_gr_arb_vec_dot(arb_t res, const arb_t initial, int subtract, arb_srcptr vec1, arb_srcptr vec2, slong len, gr_ctx_t ctx)
{
    arb_dot(res, initial, subtract, vec1, 1, vec2, 1, len, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_vec_dot_rev(arb_t res, const arb_t initial, int subtract, arb_srcptr vec1, arb_srcptr vec2, slong len, gr_ctx_t ctx)
{
    arb_dot(res, initial, subtract, vec1, 1, vec2 + len - 1, -1, len, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_poly_mullow(arb_ptr res,
    arb_srcptr poly1, slong len1,
    arb_srcptr poly2, slong len2, slong n, gr_ctx_t ctx)
{
    _arb_poly_mullow(res, poly1, len1, poly2, len2, n, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

/* hidden feature: also works with arb ctx */
int
_gr_acb_poly_roots(gr_vec_t roots, gr_vec_t mult, const gr_poly_t poly, int flags, gr_ctx_t ctx);

int
_gr_arb_poly_roots(gr_vec_t roots, gr_vec_t mult, const gr_poly_t poly, int flags, gr_ctx_t ctx)
{
    int status;
    slong i;
    acb_poly_t tmp;
    acb_poly_init(tmp);
    acb_poly_fit_length(tmp, poly->length);
    for (i = 0; i < poly->length; i++)
        acb_set_arb(tmp->coeffs + i, ((arb_srcptr) poly->coeffs) + i);
    _acb_poly_set_length(tmp, poly->length);
    status = _gr_acb_poly_roots(roots, mult, (gr_poly_struct *) tmp, flags, ctx);
    acb_poly_clear(tmp);
    return status;
}

int
_gr_arb_poly_roots_other(gr_vec_t roots, gr_vec_t mult, const gr_poly_t poly, gr_ctx_t other_ctx, int flags, gr_ctx_t ctx)
{
    if (poly->length == 0)
        return GR_DOMAIN;

    if (other_ctx->which_ring == GR_CTX_RR_ARB)
    {
        return _gr_arb_poly_roots(roots, mult, poly, flags, ctx);
    }

    if (other_ctx->which_ring == GR_CTX_FMPZ)
    {
        gr_ctx_t ZZ;
        slong i, j, deg, deg2;
        acb_ptr croots;
        int status = GR_SUCCESS;

        deg = poly->length - 1;

        gr_ctx_init_fmpz(ZZ);

        gr_vec_set_length(roots, 0, ctx);
        gr_vec_set_length(mult, 0, ZZ);

        if (deg != 0)
        {
            fmpz_poly_factor_t fac;
            fmpz_poly_factor_init(fac);
            fmpz_poly_factor_squarefree(fac, (const fmpz_poly_struct *) poly);

            for (i = 0; i < fac->num; i++)
            {
                deg2 = fmpz_poly_degree(fac->p + i);

                croots = _acb_vec_init(deg2);
                arb_fmpz_poly_complex_roots(croots, fac->p + i, 0, ARB_CTX_PREC(ctx));

                for (j = 0; j < deg2; j++)
                {
                    if (acb_is_real(croots + j))
                    {
                        fmpz m2 = fac->exp[i];
                        GR_MUST_SUCCEED(gr_vec_append(roots, acb_realref(croots + j), ctx));
                        GR_MUST_SUCCEED(gr_vec_append(mult, &m2, ZZ));
                    }
                }

                _acb_vec_clear(croots, deg2);
            }

            fmpz_poly_factor_clear(fac);
        }

        gr_ctx_clear(ZZ);

        return status;
    }

    return GR_UNABLE;
}

int
_gr_arb_mat_mul(arb_mat_t res, const arb_mat_t x, const arb_mat_t y, gr_ctx_t ctx)
{
    arb_mat_mul(res, x, y, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_mat_det(arb_t res, const arb_mat_t x, gr_ctx_t ctx)
{
    arb_mat_det(res, x, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int _arb_methods_initialized = 0;

gr_static_method_table _arb_methods;

gr_method_tab_input _arb_methods_input[] =
{
    {GR_METHOD_CTX_WRITE,       (gr_funcptr) _gr_arb_ctx_write},
    {GR_METHOD_CTX_IS_RING,     (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_COMMUTATIVE_RING, (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_INTEGRAL_DOMAIN,  (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_FIELD,            (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_UNIQUE_FACTORIZATION_DOMAIN,
                                (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_FINITE,
                                (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_FINITE_CHARACTERISTIC,
                                (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_ALGEBRAICALLY_CLOSED,
                                (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_ORDERED_RING,
                                (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_EXACT,    (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_CANONICAL,
                                (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_HAS_REAL_PREC, (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_SET_REAL_PREC, (gr_funcptr) _gr_arb_ctx_set_real_prec},
    {GR_METHOD_CTX_GET_REAL_PREC, (gr_funcptr) _gr_arb_ctx_get_real_prec},

    {GR_METHOD_INIT,            (gr_funcptr) _gr_arb_init},
    {GR_METHOD_CLEAR,           (gr_funcptr) _gr_arb_clear},
    {GR_METHOD_SWAP,            (gr_funcptr) _gr_arb_swap},
    {GR_METHOD_RANDTEST,        (gr_funcptr) _gr_arb_randtest},
    {GR_METHOD_WRITE,           (gr_funcptr) _gr_arb_write},
    {GR_METHOD_ZERO,            (gr_funcptr) _gr_arb_zero},
    {GR_METHOD_ONE,             (gr_funcptr) _gr_arb_one},
    {GR_METHOD_IS_ZERO,         (gr_funcptr) _gr_arb_is_zero},
    {GR_METHOD_IS_ONE,          (gr_funcptr) _gr_arb_is_one},
    {GR_METHOD_IS_NEG_ONE,      (gr_funcptr) _gr_arb_is_neg_one},
    {GR_METHOD_EQUAL,           (gr_funcptr) _gr_arb_equal},
    {GR_METHOD_SET,             (gr_funcptr) _gr_arb_set},
    {GR_METHOD_SET_SI,          (gr_funcptr) _gr_arb_set_si},
    {GR_METHOD_SET_UI,          (gr_funcptr) _gr_arb_set_ui},
    {GR_METHOD_SET_FMPZ,        (gr_funcptr) _gr_arb_set_fmpz},
    {GR_METHOD_SET_FMPQ,        (gr_funcptr) _gr_arb_set_fmpq},
    {GR_METHOD_SET_STR,         (gr_funcptr) _gr_arb_set_str},
    {GR_METHOD_SET_D,           (gr_funcptr) _gr_arb_set_d},
    {GR_METHOD_SET_OTHER,       (gr_funcptr) _gr_arb_set_other},
    {GR_METHOD_GET_SI,          (gr_funcptr) _gr_arb_get_si},
    {GR_METHOD_GET_UI,          (gr_funcptr) _gr_arb_get_ui},
    {GR_METHOD_GET_FMPZ,        (gr_funcptr) _gr_arb_get_fmpz},
    {GR_METHOD_GET_D,           (gr_funcptr) _gr_arb_get_d},
    {GR_METHOD_NEG,             (gr_funcptr) _gr_arb_neg},
    {GR_METHOD_ADD,             (gr_funcptr) _gr_arb_add},
    {GR_METHOD_ADD_UI,          (gr_funcptr) _gr_arb_add_ui},
    {GR_METHOD_ADD_SI,          (gr_funcptr) _gr_arb_add_si},
    {GR_METHOD_ADD_FMPZ,        (gr_funcptr) _gr_arb_add_fmpz},
    {GR_METHOD_SUB,             (gr_funcptr) _gr_arb_sub},
    {GR_METHOD_SUB_UI,          (gr_funcptr) _gr_arb_sub_ui},
    {GR_METHOD_SUB_SI,          (gr_funcptr) _gr_arb_sub_si},
    {GR_METHOD_SUB_FMPZ,        (gr_funcptr) _gr_arb_sub_fmpz},
    {GR_METHOD_MUL,             (gr_funcptr) _gr_arb_mul},
    {GR_METHOD_MUL_UI,          (gr_funcptr) _gr_arb_mul_ui},
    {GR_METHOD_MUL_SI,          (gr_funcptr) _gr_arb_mul_si},
    {GR_METHOD_MUL_FMPZ,        (gr_funcptr) _gr_arb_mul_fmpz},
    {GR_METHOD_MUL_TWO,         (gr_funcptr) _gr_arb_mul_two},
    {GR_METHOD_ADDMUL,          (gr_funcptr) _gr_arb_addmul},
    {GR_METHOD_SUBMUL,          (gr_funcptr) _gr_arb_submul},
    {GR_METHOD_SQR,             (gr_funcptr) _gr_arb_sqr},
    {GR_METHOD_MUL_2EXP_SI,        (gr_funcptr) _gr_arb_mul_2exp_si},
    {GR_METHOD_MUL_2EXP_FMPZ,      (gr_funcptr) _gr_arb_mul_2exp_fmpz},
    {GR_METHOD_SET_FMPZ_2EXP_FMPZ, (gr_funcptr) _gr_arb_set_fmpz_2exp_fmpz},
    {GR_METHOD_GET_FMPZ_2EXP_FMPZ, (gr_funcptr) _gr_arb_get_fmpz_2exp_fmpz},
    {GR_METHOD_DIV,             (gr_funcptr) _gr_arb_div},
    {GR_METHOD_DIV_UI,          (gr_funcptr) _gr_arb_div_ui},
    {GR_METHOD_DIV_SI,          (gr_funcptr) _gr_arb_div_si},
    {GR_METHOD_DIV_FMPZ,        (gr_funcptr) _gr_arb_div_fmpz},
    {GR_METHOD_INV,             (gr_funcptr) _gr_arb_inv},
    {GR_METHOD_IS_INVERTIBLE,   (gr_funcptr) _gr_arb_is_invertible},
    {GR_METHOD_POW,             (gr_funcptr) _gr_arb_pow},
    {GR_METHOD_POW_UI,          (gr_funcptr) _gr_arb_pow_ui},
    {GR_METHOD_POW_SI,          (gr_funcptr) _gr_arb_pow_si},
    {GR_METHOD_POW_FMPZ,        (gr_funcptr) _gr_arb_pow_fmpz},
    {GR_METHOD_POW_FMPQ,        (gr_funcptr) _gr_arb_pow_fmpq},
    {GR_METHOD_IS_SQUARE,       (gr_funcptr) _gr_arb_is_square},
    {GR_METHOD_SQRT,            (gr_funcptr) _gr_arb_sqrt},
    {GR_METHOD_RSQRT,           (gr_funcptr) _gr_arb_rsqrt},
    {GR_METHOD_FLOOR,           (gr_funcptr) _gr_arb_floor},
    {GR_METHOD_CEIL,            (gr_funcptr) _gr_arb_ceil},
    {GR_METHOD_TRUNC,           (gr_funcptr) _gr_arb_trunc},
    {GR_METHOD_NINT,            (gr_funcptr) _gr_arb_nint},
    {GR_METHOD_ABS,             (gr_funcptr) _gr_arb_abs},
    {GR_METHOD_CONJ,            (gr_funcptr) _gr_arb_conj},
    {GR_METHOD_RE,              (gr_funcptr) _gr_arb_set},
    {GR_METHOD_IM,              (gr_funcptr) _gr_arb_im},
    {GR_METHOD_SGN,             (gr_funcptr) _gr_arb_sgn},
    {GR_METHOD_CSGN,            (gr_funcptr) _gr_arb_sgn},
    {GR_METHOD_CMP,             (gr_funcptr) _gr_arb_cmp},
    {GR_METHOD_CMPABS,          (gr_funcptr) _gr_arb_cmpabs},
    {GR_METHOD_I,               (gr_funcptr) gr_not_in_domain},
    {GR_METHOD_PI,              (gr_funcptr) _gr_arb_pi},
    {GR_METHOD_EULER,           (gr_funcptr) _gr_arb_euler},
    {GR_METHOD_CATALAN,         (gr_funcptr) _gr_arb_catalan},
    {GR_METHOD_KHINCHIN,        (gr_funcptr) _gr_arb_khinchin},
    {GR_METHOD_GLAISHER,        (gr_funcptr) _gr_arb_glaisher},
    {GR_METHOD_EXP,             (gr_funcptr) _gr_arb_exp},
    {GR_METHOD_LOG,             (gr_funcptr) _gr_arb_log},
    {GR_METHOD_SIN,             (gr_funcptr) _gr_arb_sin},
    {GR_METHOD_COS,             (gr_funcptr) _gr_arb_cos},
    {GR_METHOD_TAN,             (gr_funcptr) _gr_arb_tan},
    {GR_METHOD_ATAN,            (gr_funcptr) _gr_arb_atan},
    {GR_METHOD_BERNOULLI_UI,    (gr_funcptr) _gr_arb_bernoulli_ui},
    {GR_METHOD_BERNOULLI_FMPZ,  (gr_funcptr) _gr_arb_bernoulli_fmpz},
    {GR_METHOD_ERF,             (gr_funcptr) _gr_arb_erf},
    {GR_METHOD_ERFI,            (gr_funcptr) _gr_arb_erfi},
    {GR_METHOD_ERFC,            (gr_funcptr) _gr_arb_erfc},
    {GR_METHOD_GAMMA,           (gr_funcptr) _gr_arb_gamma},
    {GR_METHOD_RGAMMA,          (gr_funcptr) _gr_arb_rgamma},
    {GR_METHOD_LGAMMA,          (gr_funcptr) _gr_arb_lgamma},
    {GR_METHOD_DIGAMMA,         (gr_funcptr) _gr_arb_digamma},
    {GR_METHOD_ZETA,            (gr_funcptr) _gr_arb_zeta},
    {GR_METHOD_VEC_DOT,         (gr_funcptr) _gr_arb_vec_dot},
    {GR_METHOD_VEC_DOT_REV,     (gr_funcptr) _gr_arb_vec_dot_rev},
    {GR_METHOD_POLY_MULLOW,     (gr_funcptr) _gr_arb_poly_mullow},
    {GR_METHOD_POLY_ROOTS,      (gr_funcptr) _gr_arb_poly_roots},
    {GR_METHOD_POLY_ROOTS_OTHER,(gr_funcptr) _gr_arb_poly_roots_other},
    {GR_METHOD_MAT_MUL,         (gr_funcptr) _gr_arb_mat_mul},
    {GR_METHOD_MAT_DET,         (gr_funcptr) _gr_arb_mat_det},
    {0,                         (gr_funcptr) NULL},
};

void
gr_ctx_init_real_arb(gr_ctx_t ctx, slong prec)
{
    ctx->which_ring = GR_CTX_RR_ARB;
    ctx->sizeof_elem = sizeof(arb_struct);
    ctx->size_limit = WORD_MAX;

    ARB_CTX_PREC(ctx) = FLINT_MAX(2, FLINT_MIN(prec, WORD_MAX / 8));
    ctx->methods = _arb_methods;

    if (!_arb_methods_initialized)
    {
        gr_method_tab_init(_arb_methods, _arb_methods_input);
        _arb_methods_initialized = 1;
    }
}
