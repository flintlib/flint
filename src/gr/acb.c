/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly_factor.h"
#include "acb.h"
#include "acb_poly.h"
#include "acb_mat.h"
#include "arb_hypgeom.h"
#include "acb_hypgeom.h"
#include "acb_dirichlet.h"
#include "acb_modular.h"
#include "acb_elliptic.h"
#include "acf.h"
#include "fmpzi.h"
#include "qqbar.h"
#include "arb_fmpz_poly.h"
#include "gr.h"
#include "gr_generic.h"
#include "gr_vec.h"
#include "gr_poly.h"

typedef struct
{
    slong prec;
}
gr_acb_ctx;

#define ACB_CTX_PREC(ring_ctx) (((gr_acb_ctx *)((ring_ctx)))->prec)

int _gr_acb_ctx_set_real_prec(gr_ctx_t ctx, slong prec)
{
    prec = FLINT_MAX(prec, 2);
    prec = FLINT_MIN(prec, WORD_MAX / 8);

    ACB_CTX_PREC(ctx) = prec;
    return GR_SUCCESS;
}

int _gr_acb_ctx_get_real_prec(slong * res, gr_ctx_t ctx)
{
    *res = ACB_CTX_PREC(ctx);
    return GR_SUCCESS;
}

int
_gr_acb_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    gr_stream_write(out, "Complex numbers (acb, prec = ");
    gr_stream_write_si(out, ACB_CTX_PREC(ctx));
    gr_stream_write(out, ")");
    return GR_SUCCESS;
}

void
_gr_acb_init(acb_t x, const gr_ctx_t ctx)
{
    acb_init(x);
}

void
_gr_acb_clear(acb_t x, const gr_ctx_t ctx)
{
    acb_clear(x);
}

void
_gr_acb_swap(acb_t x, acb_t y, const gr_ctx_t ctx)
{
    acb_t t;
    *t = *x;
    *x = *y;
    *y = *t;
}

void
_gr_acb_set_shallow(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    *res = *x;
}

/* todo: limits */
int
_gr_acb_randtest(acb_t res, flint_rand_t state, const gr_ctx_t ctx)
{
    acb_randtest(res, state, ACB_CTX_PREC(ctx), 10);
    return GR_SUCCESS;
}

int
__gr_acb_write(gr_stream_t out, const acb_t x, slong digits, int flags, const gr_ctx_t ctx)
{
    if (arb_is_zero(acb_imagref(x)))
    {
        gr_stream_write_free(out, arb_get_str(acb_realref(x), digits, flags));
    }
    else if (arb_is_zero(acb_realref(x)))
    {
        gr_stream_write_free(out, arb_get_str(acb_imagref(x), digits, flags));
        gr_stream_write(out, "*I");
    }
    else
    {
        gr_stream_write(out, "(");
        gr_stream_write_free(out, arb_get_str(acb_realref(x), digits, flags));

        if ((arb_is_exact(acb_imagref(x)) || (flags & ARB_STR_NO_RADIUS))
                && arf_sgn(arb_midref(acb_imagref(x))) < 0)
        {
            arb_t t;
            arb_init(t);
            arb_neg(t, acb_imagref(x));
            gr_stream_write(out, " - ");
            gr_stream_write_free(out, arb_get_str(t, digits, flags));
            arb_clear(t);
        }
        else
        {
            gr_stream_write(out, " + ");
            gr_stream_write_free(out, arb_get_str(acb_imagref(x), digits, flags));
        }

        gr_stream_write(out, "*I)");
    }

    return GR_SUCCESS;
}
int
_gr_acb_write(gr_stream_t out, const acb_t x, const gr_ctx_t ctx)
{
    return __gr_acb_write(out, x, ACB_CTX_PREC(ctx) * 0.30102999566398 + 1, 0, ctx);
}

int
_gr_acb_write_n(gr_stream_t out, gr_srcptr x, slong n, gr_ctx_t ctx)
{
    return __gr_acb_write(out, x, FLINT_MAX(n, 1), ARB_STR_NO_RADIUS, ctx);
}

int
_gr_acb_zero(acb_t x, const gr_ctx_t ctx)
{
    acb_zero(x);
    return GR_SUCCESS;
}

int
_gr_acb_one(acb_t x, const gr_ctx_t ctx)
{
    acb_one(x);
    return GR_SUCCESS;
}

int
_gr_acb_set_si(acb_t res, slong v, const gr_ctx_t ctx)
{
    acb_set_si(res, v);
    acb_set_round(res, res, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_set_ui(acb_t res, ulong v, const gr_ctx_t ctx)
{
    acb_set_ui(res, v);
    acb_set_round(res, res, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_set_fmpz(acb_t res, const fmpz_t v, const gr_ctx_t ctx)
{
    acb_set_round_fmpz(res, v, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_set_fmpq(acb_t res, const fmpq_t v, const gr_ctx_t ctx)
{
    acb_set_fmpq(res, v, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_set_d(acb_t res, double x, const gr_ctx_t ctx)
{
    acb_set_d(res, x);
    arb_set_round(acb_realref(res), acb_realref(res), ACB_CTX_PREC(ctx));

    if (!arb_is_finite(acb_realref(res)))
        return GR_DOMAIN;

    return GR_SUCCESS;
}

int
_gr_ca_get_acb_with_prec(acb_t res, gr_srcptr x, gr_ctx_t x_ctx, slong prec);

int
_gr_acb_set_other(acb_t res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
{
    switch (x_ctx->which_ring)
    {
        case GR_CTX_FMPZ:
            return _gr_acb_set_fmpz(res, x, ctx);

        case GR_CTX_FMPQ:
            return _gr_acb_set_fmpq(res, x, ctx);

        case GR_CTX_FMPZI:
            arb_set_round_fmpz(acb_realref(res), fmpzi_realref((const fmpzi_struct *) x), ACB_CTX_PREC(ctx));
            arb_set_round_fmpz(acb_imagref(res), fmpzi_imagref((const fmpzi_struct *) x), ACB_CTX_PREC(ctx));
            return GR_SUCCESS;

        case GR_CTX_REAL_ALGEBRAIC_QQBAR:
        case GR_CTX_COMPLEX_ALGEBRAIC_QQBAR:
            qqbar_get_acb(res, x, ACB_CTX_PREC(ctx));
            return GR_SUCCESS;

        case GR_CTX_RR_CA:
        case GR_CTX_REAL_ALGEBRAIC_CA:
        case GR_CTX_CC_CA:
        case GR_CTX_COMPLEX_ALGEBRAIC_CA:
            return _gr_ca_get_acb_with_prec(res, x, x_ctx, ACB_CTX_PREC(ctx));

        case GR_CTX_REAL_FLOAT_ARF:
            if (arf_is_finite(x))
            {
                arb_set_arf(acb_realref(res), x);
                arb_set_round(acb_realref(res), acb_realref(res), ACB_CTX_PREC(ctx));
                arb_zero(acb_imagref(res));
                return GR_SUCCESS;
            }
            else
            {
                return GR_DOMAIN;
            }

        case GR_CTX_COMPLEX_FLOAT_ACF:
            if (arf_is_finite(acf_realref((acf_srcptr) x)) && arf_is_finite(acf_imagref((acf_srcptr) x)))
            {
                arb_set_arf(acb_realref(res), acf_realref((acf_srcptr) x));
                arb_set_arf(acb_imagref(res), acf_imagref((acf_srcptr) x));
                acb_set_round(res, res, ACB_CTX_PREC(ctx));
                return GR_SUCCESS;
            }
            else
            {
                return GR_DOMAIN;
            }

        case GR_CTX_RR_ARB:
            arb_set_round(acb_realref(res), x, ACB_CTX_PREC(ctx));
            arb_zero(acb_imagref(res));
            return GR_SUCCESS;

        case GR_CTX_CC_ACB:
            acb_set_round(res, x, ACB_CTX_PREC(ctx));
            return GR_SUCCESS;
    }

    return gr_generic_set_other(res, x, x_ctx, ctx);
}

/* xxx: assumes that ctx are not read */
int _gr_arf_get_fmpz(fmpz_t res, const arf_t x, const gr_ctx_t ctx);
int _gr_arf_get_si(slong * res, const arf_t x, const gr_ctx_t ctx);
int _gr_arf_get_ui(ulong * res, const arf_t x, const gr_ctx_t ctx);

int
_gr_acb_get_fmpz(fmpz_t res, const acb_t x, const gr_ctx_t ctx)
{
    if (!acb_is_int(x))
    {
        if (acb_contains_int(x))
            return GR_UNABLE;
        else
            return GR_DOMAIN;
    }

    return _gr_arf_get_fmpz(res, arb_midref(acb_realref(x)), NULL);
}

int
_gr_acb_get_si(slong * res, const acb_t x, const gr_ctx_t ctx)
{
    if (!acb_is_int(x))
    {
        if (acb_contains_int(x))
            return GR_UNABLE;
        else
            return GR_DOMAIN;
    }

    return _gr_arf_get_si(res, arb_midref(acb_realref(x)), NULL);
}

int
_gr_acb_get_ui(ulong * res, const acb_t x, const gr_ctx_t ctx)
{
    if (!acb_is_int(x))
    {
        if (acb_contains_int(x))
            return GR_UNABLE;
        else
            return GR_DOMAIN;
    }

    return _gr_arf_get_ui(res, arb_midref(acb_realref(x)), NULL);
}

int
_gr_acb_get_d(double * res, const acb_t x, const gr_ctx_t ctx)
{
    *res = arf_get_d(arb_midref(acb_realref(x)), ARF_RND_NEAR);
    return GR_SUCCESS;
}

truth_t
_gr_acb_is_zero(const acb_t x, const gr_ctx_t ctx)
{
    if (acb_is_zero(x))
        return T_TRUE;

    if (acb_contains_zero(x))
        return T_UNKNOWN;

    return T_FALSE;
}

truth_t
_gr_acb_is_one(const acb_t x, const gr_ctx_t ctx)
{
    if (acb_is_one(x))
        return T_TRUE;

    if (arb_contains_zero(acb_imagref(x)) && arb_contains_si(acb_realref(x), 1))
        return T_UNKNOWN;

    return T_FALSE;
}

truth_t
_gr_acb_is_neg_one(const acb_t x, const gr_ctx_t ctx)
{
    if (acb_equal_si(x, -1))
        return T_TRUE;

    if (arb_contains_zero(acb_imagref(x)) && arb_contains_si(acb_realref(x), -1))
        return T_UNKNOWN;

    return T_FALSE;
}

truth_t
_gr_acb_equal(const acb_t x, const acb_t y, const gr_ctx_t ctx)
{
    if (acb_is_exact(x) && acb_equal(x, y))
        return T_TRUE;

    if (acb_overlaps(x, y))
        return T_UNKNOWN;

    return T_FALSE;
}

int
_gr_acb_set(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    acb_set(res, x);
    return GR_SUCCESS;
}

int
_gr_acb_neg(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    acb_neg(res, x);
    return GR_SUCCESS;
}

int
_gr_acb_add(acb_t res, const acb_t x, const acb_t y, const gr_ctx_t ctx)
{
    acb_add(res, x, y, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_add_si(acb_t res, const acb_t x, slong y, const gr_ctx_t ctx)
{
    acb_add_si(res, x, y, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_add_ui(acb_t res, const acb_t x, ulong y, const gr_ctx_t ctx)
{
    acb_add_ui(res, x, y, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_add_fmpz(acb_t res, const acb_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    acb_add_fmpz(res, x, y, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_sub(acb_t res, const acb_t x, const acb_t y, const gr_ctx_t ctx)
{
    acb_sub(res, x, y, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_sub_si(acb_t res, const acb_t x, slong y, const gr_ctx_t ctx)
{
    acb_sub_si(res, x, y, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_sub_ui(acb_t res, const acb_t x, ulong y, const gr_ctx_t ctx)
{
    acb_sub_ui(res, x, y, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_sub_fmpz(acb_t res, const acb_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    acb_sub_fmpz(res, x, y, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_mul(acb_t res, const acb_t x, const acb_t y, const gr_ctx_t ctx)
{
    acb_mul(res, x, y, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_mul_si(acb_t res, const acb_t x, slong y, const gr_ctx_t ctx)
{
    acb_mul_si(res, x, y, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_mul_ui(acb_t res, const acb_t x, ulong y, const gr_ctx_t ctx)
{
    acb_mul_ui(res, x, y, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_mul_fmpz(acb_t res, const acb_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    acb_mul_fmpz(res, x, y, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_addmul(acb_t res, const acb_t x, const acb_t y, const gr_ctx_t ctx)
{
    acb_addmul(res, x, y, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_submul(acb_t res, const acb_t x, const acb_t y, const gr_ctx_t ctx)
{
    acb_submul(res, x, y, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_mul_two(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    acb_mul_2exp_si(res, x, 1);
    return GR_SUCCESS;
}

int
_gr_acb_sqr(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    acb_sqr(res, x, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_mul_2exp_si(acb_t res, const acb_t x, slong y, const gr_ctx_t ctx)
{
    acb_mul_2exp_si(res, x, y);
    return GR_SUCCESS;
}

int
_gr_acb_mul_2exp_fmpz(acb_t res, const acb_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    acb_mul_2exp_fmpz(res, x, y);
    return GR_SUCCESS;
}

int
_gr_acb_set_fmpz_2exp_fmpz(acb_t res, const fmpz_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    arb_set_fmpz_2exp(acb_realref(res), x, y);
    arb_zero(acb_imagref(res));
    return GR_SUCCESS;
}

int
_gr_acb_get_fmpz_2exp_fmpz(fmpz_t res1, fmpz_t res2, const acb_t x, const gr_ctx_t ctx)
{
    if (!acb_is_exact(x) || !acb_is_real(x))
        return GR_UNABLE;

    if (!acb_is_finite(x))
        return GR_DOMAIN;

    arf_get_fmpz_2exp(res1, res2, arb_midref(acb_realref(x)));
    return GR_SUCCESS;
}

int
_gr_acb_inv(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    if (acb_is_zero(x))
    {
        return GR_DOMAIN;
    }
    else
    {
        acb_inv(res, x, ACB_CTX_PREC(ctx));
        if (acb_is_finite(res))
            return GR_SUCCESS;
        else
            return GR_UNABLE;
    }
}

int
_gr_acb_div(acb_t res, const acb_t x, const acb_t y, const gr_ctx_t ctx)
{
    if (acb_is_zero(y))
    {
        return GR_DOMAIN;
    }
    else
    {
        acb_div(res, x, y, ACB_CTX_PREC(ctx));

        if (acb_is_finite(res))
            return GR_SUCCESS;
        else
            return GR_UNABLE;
    }
}

int
_gr_acb_div_si(acb_t res, const acb_t x, slong y, const gr_ctx_t ctx)
{
    if (y == 0)
    {
        return GR_DOMAIN;
    }
    else
    {
        acb_div_si(res, x, y, ACB_CTX_PREC(ctx));
        return GR_SUCCESS;
    }
}

int
_gr_acb_div_ui(acb_t res, const acb_t x, ulong y, const gr_ctx_t ctx)
{
    if (y == 0)
    {
        return GR_DOMAIN;
    }
    else
    {
        acb_div_ui(res, x, y, ACB_CTX_PREC(ctx));
        return GR_SUCCESS;
    }
}

int
_gr_acb_div_fmpz(acb_t res, const acb_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    if (fmpz_is_zero(y))
    {
        return GR_DOMAIN;
    }
    else
    {
        acb_div_fmpz(res, x, y, ACB_CTX_PREC(ctx));
        return GR_SUCCESS;
    }
}

truth_t
_gr_acb_is_invertible(const acb_t x, const gr_ctx_t ctx)
{
    if (acb_is_zero(x))
        return T_FALSE;

    if (acb_contains_zero(x))
        return T_UNKNOWN;

    return T_TRUE;
}

int
_gr_acb_pow_ui(acb_t res, const acb_t x, ulong exp, const gr_ctx_t ctx)
{
    acb_pow_ui(res, x, exp, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_pow_si(acb_t res, const acb_t x, slong exp, const gr_ctx_t ctx)
{
    if (exp < 0 && acb_is_zero(x))
    {
        return GR_DOMAIN;
    }
    else if (exp < 0 && acb_contains_zero(x))
    {
        return GR_UNABLE;
    }
    else
    {
        fmpz_t t;
        fmpz_init_set_si(t, exp);
        acb_pow_fmpz(res, x, t, ACB_CTX_PREC(ctx));
        fmpz_clear(t);
        return GR_SUCCESS;
    }
}

int
_gr_acb_pow_fmpz(acb_t res, const acb_t x, const fmpz_t exp, const gr_ctx_t ctx)
{
    if (fmpz_sgn(exp) < 0 && acb_is_zero(x))
    {
        return GR_DOMAIN;
    }
    else if (fmpz_sgn(exp) < 0 && acb_contains_zero(x))
    {
        return GR_UNABLE;
    }
    else
    {
        acb_pow_fmpz(res, x, exp, ACB_CTX_PREC(ctx));
        return GR_SUCCESS;
    }
}

int
_gr_acb_pow(acb_t res, const acb_t x, const acb_t exp, const gr_ctx_t ctx)
{
    if (acb_is_int(exp))
    {
        if (arf_sgn(arb_midref(acb_realref(exp))) < 0)
        {
            if (acb_is_zero(x))
                return GR_DOMAIN;

            if (acb_contains_zero(x))
                return GR_UNABLE;
        }

        acb_pow(res, x, exp, ACB_CTX_PREC(ctx));
        return GR_SUCCESS;
    }
    else if (!acb_contains_zero(x) || arb_is_positive(acb_realref(exp)))
    {
        acb_pow(res, x, exp, ACB_CTX_PREC(ctx));

        if (!acb_is_finite(res))
            return GR_UNABLE;

        return GR_SUCCESS;
    }
    else if (acb_is_zero(x) && arb_is_negative(acb_realref(exp)))
    {
        return GR_DOMAIN;
    }
    else
    {
        return GR_UNABLE;
    }
}

int
_gr_acb_pow_fmpq(acb_t res, const acb_t x, const fmpq_t exp, const gr_ctx_t ctx)
{
    int status;
    acb_t t;
    acb_init(t);
    acb_set_fmpq(t, exp, ACB_CTX_PREC(ctx) + 20);
    status = _gr_acb_pow(res, x, t, ctx);
    acb_clear(t);
    return status;
}

truth_t
_gr_acb_is_square(const acb_t x, const gr_ctx_t ctx)
{
    return T_TRUE;
}

int
_gr_acb_sqrt(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    acb_sqrt(res, x, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_rsqrt(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    if (!acb_contains_zero(x))
    {
        acb_rsqrt(res, x, ACB_CTX_PREC(ctx));
        return GR_SUCCESS;
    }
    else if (acb_is_zero(x))
    {
        return GR_DOMAIN;
    }
    else
    {
        return GR_UNABLE;
    }
}

int
_gr_acb_floor(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    arb_floor(acb_realref(res), acb_realref(x), ACB_CTX_PREC(ctx));
    arb_zero(acb_imagref(res));
    return GR_SUCCESS;
}

int
_gr_acb_ceil(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    arb_ceil(acb_realref(res), acb_realref(x), ACB_CTX_PREC(ctx));
    arb_zero(acb_imagref(res));
    return GR_SUCCESS;
}

int
_gr_acb_trunc(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    arb_trunc(acb_realref(res), acb_realref(x), ACB_CTX_PREC(ctx));
    arb_zero(acb_imagref(res));
    return GR_SUCCESS;
}

int
_gr_acb_nint(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    arb_zero(acb_imagref(res));
    arb_nint(acb_realref(res), acb_realref(x), ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_i(acb_t res, const gr_ctx_t ctx)
{
    acb_onei(res);
    return GR_SUCCESS;
}

int
_gr_acb_abs(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    acb_abs(acb_realref(res), x, ACB_CTX_PREC(ctx));
    arb_zero(acb_imagref(res));
    return GR_SUCCESS;
}

int
_gr_acb_conj(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    acb_conj(res, x);
    return GR_SUCCESS;
}

int
_gr_acb_re(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    arb_set(acb_realref(res), acb_realref(x));
    arb_zero(acb_imagref(res));
    return GR_SUCCESS;
}

int
_gr_acb_im(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    arb_set(acb_realref(res), acb_imagref(x));
    arb_zero(acb_imagref(res));
    return GR_SUCCESS;
}

#define DEF_FUNC(fname) \
int \
_gr_acb_ ## fname(acb_t res, const acb_t x, const gr_ctx_t ctx) \
{ \
    acb_ ## fname(res, x, ACB_CTX_PREC(ctx)); \
    return GR_SUCCESS; \
} \

#define DEF_2FUNC(fname) \
int \
_gr_acb_ ## fname(acb_t res1, acb_t res2, const acb_t x, const gr_ctx_t ctx) \
{ \
    acb_ ## fname(res1, res2, x, ACB_CTX_PREC(ctx)); \
    return GR_SUCCESS; \
} \

#define DEF_FUNC2(fname) \
int \
_gr_acb_ ## fname(acb_t res, const acb_t x, const acb_t y, const gr_ctx_t ctx) \
{ \
    acb_ ## fname(res, x, y, ACB_CTX_PREC(ctx)); \
    return GR_SUCCESS; \
} \

#define DEF_FUNC_SING(fname) \
int \
_gr_acb_ ## fname(acb_t res, const acb_t x, const gr_ctx_t ctx) \
{ \
    acb_ ## fname(res, x, ACB_CTX_PREC(ctx)); \
    return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE; \
} \

#define DEF_FUNC2_SING(fname) \
int \
_gr_acb_ ## fname(acb_t res, const acb_t x, const acb_t y, const gr_ctx_t ctx) \
{ \
    acb_ ## fname(res, x, y, ACB_CTX_PREC(ctx)); \
    return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE; \
} \

DEF_FUNC(sgn)

int
_gr_acb_csgn(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    acb_csgn(acb_realref(res), x);
    arb_zero(acb_imagref(res));
    return GR_SUCCESS;
}

int
_gr_acb_arg(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    acb_arg(acb_realref(res), x, ACB_CTX_PREC(ctx));
    arb_zero(acb_imagref(res));
    return GR_SUCCESS;
}

int
_gr_acb_pi(acb_t res, const gr_ctx_t ctx)
{
    acb_const_pi(res, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

DEF_FUNC(exp)
DEF_FUNC(expm1)
DEF_FUNC(exp_pi_i)

int
_gr_acb_log(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    if (acb_contains_zero(x))
    {
        if (acb_is_zero(x))
            return GR_DOMAIN;
        return GR_UNABLE;
    }

    acb_log(res, x, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_log1p(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    if (arb_contains_si(acb_realref(x), -1) && arb_contains_zero(acb_imagref(x)))
    {
        if (arb_equal_si(acb_realref(x), -1) && arb_is_zero(acb_imagref(x)))
            return GR_DOMAIN;
        return GR_UNABLE;
    }

    acb_log1p(res, x, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_log_pi_i(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    acb_t t;

    if (acb_contains_zero(x))
    {
        if (acb_is_zero(x))
            return GR_DOMAIN;
        return GR_UNABLE;
    }

    acb_log(res, x, ACB_CTX_PREC(ctx));

    acb_init(t);
    acb_const_pi(t, ACB_CTX_PREC(ctx));
    acb_mul_onei(t, t);
    acb_div(res, res, t, ACB_CTX_PREC(ctx));
    acb_clear(t);

    return GR_SUCCESS;
}


DEF_FUNC(sin)
DEF_FUNC(sin_pi)
DEF_FUNC(cos)
DEF_FUNC(cos_pi)
DEF_FUNC_SING(tan)
DEF_FUNC_SING(cot)

DEF_FUNC(sinc)
DEF_FUNC(sinc_pi)

DEF_FUNC(sinh)
DEF_FUNC(cosh)
DEF_FUNC_SING(tanh)
DEF_FUNC_SING(coth)

DEF_FUNC(asin)
DEF_FUNC(acos)

DEF_2FUNC(sin_cos)
DEF_2FUNC(sin_cos_pi)

DEF_2FUNC(sinh_cosh)


int
_gr_acb_atan(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    if (!arb_is_zero(acb_imagref(x)) && arb_contains_zero(acb_realref(x)))
    {
        if (arb_is_zero(acb_realref(x)) && (arb_is_one(acb_imagref(x)) || arb_equal_si(acb_imagref(x), -1)))
            return GR_DOMAIN;

        if (arb_contains_si(acb_imagref(x), 1) || arb_contains_si(acb_imagref(x), -1))
            return GR_UNABLE;
    }

    acb_atan(res, x, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

DEF_FUNC(asinh)
DEF_FUNC(acosh)

DEF_FUNC_SING(atanh)

int
_gr_acb_lambertw(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    fmpz_t k;
    fmpz_init(k);
    acb_lambertw(res, x, k, 0, ACB_CTX_PREC(ctx));
    fmpz_clear(k);
    return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int
_gr_acb_lambertw_fmpz(acb_t res, const acb_t x, const fmpz_t k, const gr_ctx_t ctx)
{
    acb_lambertw(res, x, k, 0, ACB_CTX_PREC(ctx));
    return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int
_gr_acb_erf(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    acb_hypgeom_erf(res, x, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_erfc(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    acb_hypgeom_erfc(res, x, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_erfi(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    acb_hypgeom_erfi(res, x, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_erfinv(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    if (acb_is_real(x))
    {
        arb_hypgeom_erfinv(acb_realref(res), acb_realref(x), ACB_CTX_PREC(ctx));
        arb_zero(acb_imagref(res));
    }
    else
    {
        acb_indeterminate(res);
    }
    return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int
_gr_acb_erfcinv(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    if (acb_is_real(x))
    {
        arb_hypgeom_erfcinv(acb_realref(res), acb_realref(x), ACB_CTX_PREC(ctx));
        arb_zero(acb_imagref(res));
    }
    else
    {
        acb_indeterminate(res);
    }
    return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int
_gr_acb_fresnel_s(acb_t res, const acb_t x, int normalized, const gr_ctx_t ctx)
{
    acb_hypgeom_fresnel(res, NULL, x, normalized, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_fresnel_c(acb_t res, const acb_t x, int normalized, const gr_ctx_t ctx)
{
    acb_hypgeom_fresnel(NULL, res, x, normalized, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_fresnel(acb_t res1, acb_t res2, const acb_t x, int normalized, const gr_ctx_t ctx)
{
    acb_hypgeom_fresnel(res1, res2, x, normalized, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_gamma_upper(acb_t res, const acb_t x, const acb_t y, int regularized, const gr_ctx_t ctx)
{
    acb_hypgeom_gamma_upper(res, x, y, regularized, ACB_CTX_PREC(ctx));
    return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int
_gr_acb_gamma_lower(acb_t res, const acb_t x, const acb_t y, int regularized, const gr_ctx_t ctx)
{
    acb_hypgeom_gamma_lower(res, x, y, regularized, ACB_CTX_PREC(ctx));
    return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int
_gr_acb_beta_lower(acb_t res, const acb_t x, const acb_t y, const acb_t z, int regularized, const gr_ctx_t ctx)
{
    acb_hypgeom_beta_lower(res, x, y, z, regularized, ACB_CTX_PREC(ctx));
    return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int
_gr_acb_exp_integral(acb_t res, const acb_t x, const acb_t y, const gr_ctx_t ctx)
{
    acb_hypgeom_expint(res, x, y, ACB_CTX_PREC(ctx));
    return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int
_gr_acb_exp_integral_ei(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    acb_hypgeom_ei(res, x, ACB_CTX_PREC(ctx));
    return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int
_gr_acb_sin_integral(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    acb_hypgeom_si(res, x, ACB_CTX_PREC(ctx));
    return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int
_gr_acb_cos_integral(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    acb_hypgeom_ci(res, x, ACB_CTX_PREC(ctx));
    return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int
_gr_acb_sinh_integral(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    acb_hypgeom_shi(res, x, ACB_CTX_PREC(ctx));
    return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int
_gr_acb_cosh_integral(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    acb_hypgeom_chi(res, x, ACB_CTX_PREC(ctx));
    return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int
_gr_acb_log_integral(acb_t res, const acb_t x, int offset, const gr_ctx_t ctx)
{
    acb_hypgeom_li(res, x, offset, ACB_CTX_PREC(ctx));
    return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int
_gr_acb_dilog(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    acb_hypgeom_dilog(res, x, ACB_CTX_PREC(ctx));
    return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int
_gr_acb_gamma(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    if (acb_is_int(x) && arb_is_nonpositive(acb_realref(x)))
    {
        return GR_DOMAIN;
    }
    else
    {
        acb_gamma(res, x, ACB_CTX_PREC(ctx));
        return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
    }
}

int
_gr_acb_gamma_fmpz(acb_t res, const fmpz_t x, const gr_ctx_t ctx)
{
    if (fmpz_sgn(x) > 0)
    {
        arb_gamma_fmpz(acb_realref(res), x, ACB_CTX_PREC(ctx));
        arb_zero(acb_imagref(res));
        return GR_SUCCESS;
    }
    else
    {
        return GR_DOMAIN;
    }
}

int
_gr_acb_gamma_fmpq(acb_t res, const fmpq_t x, const gr_ctx_t ctx)
{
    if (!fmpz_is_one(fmpq_denref(x)) || fmpz_sgn(fmpq_numref(x)) > 0)
    {
        arb_gamma_fmpq(acb_realref(res), x, ACB_CTX_PREC(ctx));
        arb_zero(acb_imagref(res));
        return GR_SUCCESS;
    }
    else
    {
        return GR_DOMAIN;
    }
}


int
_gr_acb_rgamma(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    acb_rgamma(res, x, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_lgamma(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    if (acb_is_int(x) && arb_is_nonpositive(acb_realref(x)))
    {
        return GR_DOMAIN;
    }
    else
    {
        acb_lgamma(res, x, ACB_CTX_PREC(ctx));
        return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
    }
}

int
_gr_acb_digamma(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    if (acb_is_int(x) && arb_is_nonpositive(acb_realref(x)))
    {
        return GR_DOMAIN;
    }
    else
    {
        acb_digamma(res, x, ACB_CTX_PREC(ctx));
        return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
    }
}

int
_gr_acb_fac_ui(acb_t res, ulong x, const gr_ctx_t ctx)
{
    arb_fac_ui(acb_realref(res), x, ACB_CTX_PREC(ctx));
    arb_zero(acb_imagref(res));
    return GR_SUCCESS;
}

int
_gr_acb_fac_fmpz(acb_t res, const fmpz_t x, const gr_ctx_t ctx)
{
    int status;
    fmpz_t t;
    fmpz_init(t);
    fmpz_add_ui(t, x, 1);
    status = _gr_acb_gamma_fmpz(res, t, ctx);
    fmpz_clear(t);
    return status;
}

int
_gr_acb_rising_ui(acb_t res, const acb_t x, ulong y, const gr_ctx_t ctx)
{
    acb_rising_ui(res, x, y, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_rising(acb_t res, const acb_t x, const acb_t y, const gr_ctx_t ctx)
{
    acb_rising(res, x, y, ACB_CTX_PREC(ctx));

    if (acb_is_finite(res))
        return GR_SUCCESS;
    else
        return GR_UNABLE;
}

int
_gr_acb_barnes_g(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    if (acb_is_int(x) && arb_is_nonpositive(acb_realref(x)))
    {
        return GR_DOMAIN;
    }
    else
    {
        acb_barnes_g(res, x, ACB_CTX_PREC(ctx));
        return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
    }
}

int
_gr_acb_log_barnes_g(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    if (acb_is_int(x) && arb_is_nonpositive(acb_realref(x)))
    {
        return GR_DOMAIN;
    }
    else
    {
        acb_log_barnes_g(res, x, ACB_CTX_PREC(ctx));
        return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
    }
}

int _gr_acb_bessel_j(acb_t res, const acb_t x, const acb_t y, const gr_ctx_t ctx) { acb_hypgeom_bessel_j(res, x, y, ACB_CTX_PREC(ctx)); return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE; }
int _gr_acb_bessel_y(acb_t res, const acb_t x, const acb_t y, const gr_ctx_t ctx) { acb_hypgeom_bessel_y(res, x, y, ACB_CTX_PREC(ctx)); return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE; }
int _gr_acb_bessel_i(acb_t res, const acb_t x, const acb_t y, const gr_ctx_t ctx) { acb_hypgeom_bessel_i(res, x, y, ACB_CTX_PREC(ctx)); return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE; }
int _gr_acb_bessel_k(acb_t res, const acb_t x, const acb_t y, const gr_ctx_t ctx) { acb_hypgeom_bessel_k(res, x, y, ACB_CTX_PREC(ctx)); return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE; }

int _gr_acb_bessel_j_y(acb_t res1, acb_t res2, const acb_t x, const acb_t y, const gr_ctx_t ctx)
{
    acb_hypgeom_bessel_jy(res1, res2, x, y, ACB_CTX_PREC(ctx));
    return (acb_is_finite(res1) && acb_is_finite(res2)) ? GR_SUCCESS : GR_UNABLE;
}

int _gr_acb_bessel_i_scaled(acb_t res, const acb_t x, const acb_t y, const gr_ctx_t ctx) { acb_hypgeom_bessel_i_scaled(res, x, y, ACB_CTX_PREC(ctx)); return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE; }
int _gr_acb_bessel_k_scaled(acb_t res, const acb_t x, const acb_t y, const gr_ctx_t ctx) { acb_hypgeom_bessel_k_scaled(res, x, y, ACB_CTX_PREC(ctx)); return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE; }

int _gr_acb_airy(acb_t res1, acb_t res2, acb_t res3, acb_t res4, const acb_t x, const gr_ctx_t ctx) { acb_hypgeom_airy(res1, res2, res3, res4, x, ACB_CTX_PREC(ctx)); return GR_SUCCESS; }
int _gr_acb_airy_ai(acb_t res, const acb_t x, const gr_ctx_t ctx) { acb_hypgeom_airy(res, NULL, NULL, NULL, x, ACB_CTX_PREC(ctx)); return GR_SUCCESS; }
int _gr_acb_airy_ai_prime(acb_t res, const acb_t x, const gr_ctx_t ctx) { acb_hypgeom_airy(NULL, res, NULL, NULL, x, ACB_CTX_PREC(ctx)); return GR_SUCCESS; }
int _gr_acb_airy_bi(acb_t res, const acb_t x, const gr_ctx_t ctx) { acb_hypgeom_airy(NULL, NULL, res, NULL, x, ACB_CTX_PREC(ctx)); return GR_SUCCESS; }
int _gr_acb_airy_bi_prime(acb_t res, const acb_t x, const gr_ctx_t ctx) { acb_hypgeom_airy(NULL, NULL, NULL, res, x, ACB_CTX_PREC(ctx)); return GR_SUCCESS; }

int _gr_acb_airy_ai_zero(acb_t res, const fmpz_t n, const gr_ctx_t ctx)
{
    if (fmpz_sgn(n) <= 0)
        return GR_DOMAIN;
    arb_hypgeom_airy_zero(acb_realref(res), NULL, NULL, NULL, n, ACB_CTX_PREC(ctx));
    arb_zero(acb_imagref(res));
    return GR_SUCCESS;
}

int _gr_acb_airy_bi_zero(acb_t res, const fmpz_t n, const gr_ctx_t ctx)
{
    if (fmpz_sgn(n) <= 0)
        return GR_DOMAIN;
    arb_hypgeom_airy_zero(NULL, NULL, acb_realref(res), NULL, n, ACB_CTX_PREC(ctx));
    arb_zero(acb_imagref(res));
    return GR_SUCCESS;
}

int _gr_acb_airy_ai_prime_zero(acb_t res, const fmpz_t n, const gr_ctx_t ctx)
{
    if (fmpz_sgn(n) <= 0)
        return GR_DOMAIN;
    arb_hypgeom_airy_zero(NULL, acb_realref(res), NULL, NULL, n, ACB_CTX_PREC(ctx));
    arb_zero(acb_imagref(res));
    return GR_SUCCESS;
}

int _gr_acb_airy_bi_prime_zero(acb_t res, const fmpz_t n, const gr_ctx_t ctx)
{
    if (fmpz_sgn(n) <= 0)
        return GR_DOMAIN;
    arb_hypgeom_airy_zero(NULL, NULL, NULL, acb_realref(res), n, ACB_CTX_PREC(ctx));
    arb_zero(acb_imagref(res));
    return GR_SUCCESS;
}

int _gr_acb_coulomb(acb_t res1, acb_t res2, acb_t res3, acb_t res4, const acb_t x, const acb_t y, const acb_t z, const gr_ctx_t ctx)
{
    acb_hypgeom_coulomb(res1, res2, res3, res4, x, y, z, ACB_CTX_PREC(ctx));
    return (acb_is_finite(res1) && acb_is_finite(res2) && acb_is_finite(res3) && acb_is_finite(res4)) ? GR_SUCCESS : GR_UNABLE;
}

int _gr_acb_coulomb_f(acb_t res, const acb_t x, const acb_t y, const acb_t z, const gr_ctx_t ctx) { acb_hypgeom_coulomb(res, NULL, NULL, NULL, x, y, z, ACB_CTX_PREC(ctx)); return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE; }
int _gr_acb_coulomb_g(acb_t res, const acb_t x, const acb_t y, const acb_t z, const gr_ctx_t ctx) { acb_hypgeom_coulomb(NULL, res, NULL, NULL, x, y, z, ACB_CTX_PREC(ctx)); return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE; }
int _gr_acb_coulomb_hpos(acb_t res, const acb_t x, const acb_t y, const acb_t z, const gr_ctx_t ctx) { acb_hypgeom_coulomb(NULL, NULL, res, NULL, x, y, z, ACB_CTX_PREC(ctx)); return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE; }
int _gr_acb_coulomb_hneg(acb_t res, const acb_t x, const acb_t y, const acb_t z, const gr_ctx_t ctx) { acb_hypgeom_coulomb(NULL, NULL, NULL, res, x, y, z, ACB_CTX_PREC(ctx)); return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE; }

int _gr_acb_chebyshev_t(acb_t res, const acb_t n, const acb_t x, const gr_ctx_t ctx)
{
    acb_hypgeom_chebyshev_t(res, n, x, ACB_CTX_PREC(ctx));
    return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int _gr_acb_chebyshev_u(acb_t res, const acb_t n, const acb_t x, const gr_ctx_t ctx)
{
    acb_hypgeom_chebyshev_u(res, n, x, ACB_CTX_PREC(ctx));
    return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int _gr_acb_jacobi_p(acb_t res, const acb_t n, const acb_t a, const acb_t b, const acb_t x, const gr_ctx_t ctx)
{
    acb_hypgeom_jacobi_p(res, n, a, b, x, ACB_CTX_PREC(ctx));
    return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int _gr_acb_gegenbauer_c(acb_t res, const acb_t n, const acb_t m, const acb_t x, const gr_ctx_t ctx)
{
    acb_hypgeom_gegenbauer_c(res, n, m, x, ACB_CTX_PREC(ctx));
    return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int _gr_acb_laguerre_l(acb_t res, const acb_t n, const acb_t m, const acb_t x, const gr_ctx_t ctx)
{
    acb_hypgeom_laguerre_l(res, n, m, x, ACB_CTX_PREC(ctx));
    return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int _gr_acb_hermite_h(acb_t res, const acb_t n, const acb_t x, const gr_ctx_t ctx)
{
    acb_hypgeom_hermite_h(res, n, x, ACB_CTX_PREC(ctx));
    return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int _gr_acb_legendre_p(acb_t res, const acb_t n, const acb_t m, const acb_t x, int type, const gr_ctx_t ctx)
{
    acb_hypgeom_legendre_p(res, n, m, x, type, ACB_CTX_PREC(ctx));
    return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int _gr_acb_legendre_q(acb_t res, const acb_t n, const acb_t m, const acb_t x, int type, const gr_ctx_t ctx)
{
    acb_hypgeom_legendre_q(res, n, m, x, type, ACB_CTX_PREC(ctx));
    return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int _gr_acb_legendre_p_root_ui(acb_t res, acb_t res2, ulong n, ulong k, const gr_ctx_t ctx)
{
    if (k >= n)
        return GR_DOMAIN;

    arb_hypgeom_legendre_p_ui_root(acb_realref(res), res2 ? acb_realref(res2) : NULL, n, k, ACB_CTX_PREC(ctx));
    arb_zero(acb_imagref(res));
    if (res2 != NULL)
        arb_zero(acb_imagref(res2));
    return GR_SUCCESS;
}

int _gr_acb_spherical_y_si(acb_t res, slong n, slong m, const acb_t theta, const acb_t phi, gr_ctx_t ctx)
{
    acb_hypgeom_spherical_y(res, n, m, theta, phi, ACB_CTX_PREC(ctx));
    return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int _gr_acb_hypgeom_0f1(acb_t res, const acb_t a, const acb_t x, int flags, const gr_ctx_t ctx)
{
    acb_hypgeom_0f1(res, a, x, flags, ACB_CTX_PREC(ctx));
    return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int _gr_acb_hypgeom_1f1(acb_t res, const acb_t a, const acb_t b, const acb_t x, int flags, const gr_ctx_t ctx)
{
    acb_hypgeom_1f1(res, a, b, x, flags, ACB_CTX_PREC(ctx));
    return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int _gr_acb_hypgeom_u(acb_t res, const acb_t a, const acb_t b, const acb_t x, int flags, const gr_ctx_t ctx)
{
    acb_hypgeom_u(res, a, b, x, ACB_CTX_PREC(ctx));
    return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int _gr_acb_hypgeom_2f1(acb_t res, const acb_t a, const acb_t b, const acb_t c, const acb_t x, int flags, const gr_ctx_t ctx)
{
    acb_hypgeom_2f1(res, a, b, c, x, flags, ACB_CTX_PREC(ctx));
    return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int _gr_acb_hypgeom_pfq(acb_t res, const gr_vec_t a, const gr_vec_t b, const acb_t x, int flags, const gr_ctx_t ctx)
{
    acb_hypgeom_pfq(res, a->entries, a->length, b->entries, b->length, x, flags, ACB_CTX_PREC(ctx));
    return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int
_gr_acb_zeta(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    if (arb_contains_si(acb_realref(x), 1) && arb_contains_zero(acb_imagref(x)))
    {
        if (acb_is_one(x))
            return GR_DOMAIN;
        else
            return GR_UNABLE;
    }
    else
    {
        acb_zeta(res, x, ACB_CTX_PREC(ctx));
        return GR_SUCCESS;
    }
}

int _gr_acb_hurwitz_zeta(acb_t res, const acb_t s, const acb_t a, const gr_ctx_t ctx)
{
    acb_hurwitz_zeta(res, s, a, ACB_CTX_PREC(ctx));
    return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int _gr_acb_polylog(acb_t res, const acb_t s, const acb_t z, const gr_ctx_t ctx)
{
    acb_polylog(res, s, z, ACB_CTX_PREC(ctx));
    return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int _gr_acb_polygamma(acb_t res, const acb_t s, const acb_t z, const gr_ctx_t ctx)
{
    acb_polygamma(res, s, z, ACB_CTX_PREC(ctx));
    return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int _gr_acb_lerch_phi(acb_t res, const acb_t z, const acb_t s, const acb_t a, const gr_ctx_t ctx)
{
    acb_dirichlet_lerch_phi(res, z, s, a, ACB_CTX_PREC(ctx));
    return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int _gr_acb_stieltjes(acb_t res, const fmpz_t n, const acb_t a, const gr_ctx_t ctx)
{
    acb_dirichlet_stieltjes(res, n, a, ACB_CTX_PREC(ctx));
    return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int
_gr_acb_dirichlet_eta(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    acb_dirichlet_eta(res, x, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

/* todo
int
_gr_acb_dirichlet_beta(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    acb_dirichlet_beta(res, x, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}
*/

int
_gr_acb_riemann_xi(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    acb_dirichlet_xi(res, x, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_zeta_zero(acb_t res, const fmpz_t n, const gr_ctx_t ctx)
{
    if (fmpz_sgn(n) <= 0)
    {
        return GR_UNABLE;
    }

    acb_dirichlet_zeta_zero(res, n, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_zeta_zero_vec(acb_ptr res, const fmpz_t n, slong len, const gr_ctx_t ctx)
{
    if (fmpz_sgn(n) <= 0)
    {
        return GR_UNABLE;
    }

    acb_dirichlet_zeta_zeros(res, n, len, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_zeta_nzeros(acb_t res, const acb_t t, const gr_ctx_t ctx)
{
    if (!acb_is_real(t) || !acb_is_finite(t))
        return GR_UNABLE;

    acb_dirichlet_zeta_nzeros(acb_realref(res), acb_realref(t), ACB_CTX_PREC(ctx));
    arb_zero(acb_imagref(res));
    return GR_SUCCESS;
}

int _gr_acb_modular_j(gr_ptr res, gr_srcptr tau, gr_ctx_t ctx) { acb_modular_j(res, tau, ACB_CTX_PREC(ctx)); return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE; }
int _gr_acb_modular_lambda(gr_ptr res, gr_srcptr tau, gr_ctx_t ctx) { acb_modular_lambda(res, tau, ACB_CTX_PREC(ctx)); return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE; }
int _gr_acb_modular_delta(gr_ptr res, gr_srcptr tau, gr_ctx_t ctx) { acb_modular_delta(res, tau, ACB_CTX_PREC(ctx)); return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE; }
int _gr_acb_dedekind_eta(gr_ptr res, gr_srcptr tau, gr_ctx_t ctx) { acb_modular_eta(res, tau, ACB_CTX_PREC(ctx)); return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE; }

int _gr_acb_eisenstein_g(gr_ptr res, ulong k, gr_srcptr tau, gr_ctx_t ctx)
{
    if (k == 0 || k % 2 == 1)
        return GR_DOMAIN;

    if (k == 2)
    {
        acb_t t;
        acb_init(t);
        acb_set_d(t, 0.5);
        acb_elliptic_zeta(res, t, tau, ACB_CTX_PREC(ctx));
        acb_mul_2exp_si(res, res, 1);
        acb_clear(t);
    }
    else   /* todo: better algorithm */
    {
        acb_ptr t;

        t = _acb_vec_init(k / 2 - 1);
        acb_modular_eisenstein(t, tau, k / 2 - 1, ACB_CTX_PREC(ctx));
        acb_swap((acb_ptr) res, t + (k / 2 - 2));
        _acb_vec_clear(t, k / 2 - 1);
    }

    return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int _gr_acb_eisenstein_g_vec(gr_ptr res, gr_srcptr tau, slong len, gr_ctx_t ctx)
{
    acb_modular_eisenstein(res, tau,  len, ACB_CTX_PREC(ctx));
    return _arb_vec_is_finite(res, 2 * len) ? GR_SUCCESS : GR_UNABLE;
}

int _gr_acb_eisenstein_e(gr_ptr res, ulong k, gr_srcptr tau, gr_ctx_t ctx)
{
    int status;

    status = _gr_acb_eisenstein_g(res, k, tau, ctx);

    if (status == GR_SUCCESS)
    {
        arb_t t;
        arb_init(t);
        arb_zeta_ui(t, k, ACB_CTX_PREC(ctx));
        acb_div_arb(res, res, t, ACB_CTX_PREC(ctx));
        acb_mul_2exp_si(res, res, -1);
        arb_clear(t);
    }

    return status;
}

int _gr_acb_jacobi_theta(gr_ptr res1, gr_ptr res2, gr_ptr res3, gr_ptr res4, gr_srcptr z, gr_srcptr tau, gr_ctx_t ctx)
{
    acb_modular_theta(res1, res2, res3, res4, z, tau, ACB_CTX_PREC(ctx));
    return (acb_is_finite(res1) && acb_is_finite(res2) && acb_is_finite(res3) && acb_is_finite(res4)) ? GR_SUCCESS : GR_UNABLE;
}

int _gr_acb_jacobi_theta_1(gr_ptr res, gr_srcptr z, gr_srcptr tau, gr_ctx_t ctx)
{
    acb_t t, u, v;
    acb_init(t);
    acb_init(u);
    acb_init(v);
    acb_modular_theta(res, t, u, v, z, tau, ACB_CTX_PREC(ctx));
    acb_clear(t);
    acb_clear(u);
    acb_clear(v);
    return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int _gr_acb_jacobi_theta_2(gr_ptr res, gr_srcptr z, gr_srcptr tau, gr_ctx_t ctx)
{
    acb_t t, u, v;
    acb_init(t);
    acb_init(u);
    acb_init(v);
    acb_modular_theta(t, res, u, v, z, tau, ACB_CTX_PREC(ctx));
    acb_clear(t);
    acb_clear(u);
    acb_clear(v);
    return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int _gr_acb_jacobi_theta_3(gr_ptr res, gr_srcptr z, gr_srcptr tau, gr_ctx_t ctx)
{
    acb_t t, u, v;
    acb_init(t);
    acb_init(u);
    acb_init(v);
    acb_modular_theta(t, u, res, v, z, tau, ACB_CTX_PREC(ctx));
    acb_clear(t);
    acb_clear(u);
    acb_clear(v);
    return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int _gr_acb_jacobi_theta_4(gr_ptr res, gr_srcptr z, gr_srcptr tau, gr_ctx_t ctx)
{
    acb_t t, u, v;
    acb_init(t);
    acb_init(u);
    acb_init(v);
    acb_modular_theta(t, u, v, res, z, tau, ACB_CTX_PREC(ctx));
    acb_clear(t);
    acb_clear(u);
    acb_clear(v);
    return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int _gr_acb_elliptic_k(gr_ptr res, gr_srcptr m, gr_ctx_t ctx) { acb_elliptic_k(res, m, ACB_CTX_PREC(ctx)); return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE; }
int _gr_acb_elliptic_e(gr_ptr res, gr_srcptr m, gr_ctx_t ctx) { acb_elliptic_e(res, m, ACB_CTX_PREC(ctx)); return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE; }
int _gr_acb_elliptic_pi(gr_ptr res, gr_srcptr n, gr_srcptr m, gr_ctx_t ctx) { acb_elliptic_pi(res, n, m, ACB_CTX_PREC(ctx)); return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE; }
int _gr_acb_elliptic_f(gr_ptr res, gr_srcptr phi, gr_srcptr m, int pi, gr_ctx_t ctx) { acb_elliptic_f(res, phi, m, pi, ACB_CTX_PREC(ctx)); return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE; }
int _gr_acb_elliptic_e_inc(gr_ptr res, gr_srcptr phi, gr_srcptr m, int pi, gr_ctx_t ctx) { acb_elliptic_e_inc(res, phi, m, pi, ACB_CTX_PREC(ctx)); return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE; }
int _gr_acb_elliptic_pi_inc(gr_ptr res, gr_srcptr n, gr_srcptr phi, int pi, gr_srcptr m, gr_ctx_t ctx) { acb_elliptic_pi_inc(res, n, phi, m, pi, ACB_CTX_PREC(ctx)); return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE; }

int _gr_acb_carlson_rc(gr_ptr res, gr_srcptr x, gr_srcptr y, int flags, gr_ctx_t ctx) { acb_elliptic_rf(res, x, y, y, flags, ACB_CTX_PREC(ctx)); return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE; }
int _gr_acb_carlson_rf(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_srcptr z, int flags, gr_ctx_t ctx) { acb_elliptic_rf(res, x, y, z, flags, ACB_CTX_PREC(ctx)); return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE; }
int _gr_acb_carlson_rg(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_srcptr z, int flags, gr_ctx_t ctx) { acb_elliptic_rg(res, x, y, z, flags, ACB_CTX_PREC(ctx)); return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE; }
int _gr_acb_carlson_rd(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_srcptr z, int flags, gr_ctx_t ctx) { acb_elliptic_rj(res, x, y, z, z, flags, ACB_CTX_PREC(ctx)); return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE; }
int _gr_acb_carlson_rj(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_srcptr z, gr_srcptr w, int flags, gr_ctx_t ctx) { acb_elliptic_rj(res, x, y, z, w, flags, ACB_CTX_PREC(ctx)); return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE; }

int _gr_acb_elliptic_invariants(gr_ptr res1, gr_ptr res2, gr_srcptr tau, gr_ctx_t ctx) { acb_elliptic_invariants(res1, res2, tau, ACB_CTX_PREC(ctx)); return (acb_is_finite(res1) && acb_is_finite(res2)) ? GR_SUCCESS : GR_UNABLE; }
int _gr_acb_elliptic_roots(gr_ptr res1, gr_ptr res2, gr_ptr res3, gr_srcptr tau, gr_ctx_t ctx) { acb_elliptic_roots(res1, res2, res3, tau, ACB_CTX_PREC(ctx)); return (acb_is_finite(res1) && acb_is_finite(res2) && acb_is_finite(res3)) ? GR_SUCCESS : GR_UNABLE; }

int _gr_acb_weierstrass_p(gr_ptr res, gr_srcptr z, gr_srcptr tau, gr_ctx_t ctx) { acb_elliptic_p(res, z, tau, ACB_CTX_PREC(ctx)); return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE; }
int _gr_acb_weierstrass_p_prime(gr_ptr res, gr_srcptr z, gr_srcptr tau, gr_ctx_t ctx) { acb_elliptic_p_prime(res, z, tau, ACB_CTX_PREC(ctx)); return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE; }
int _gr_acb_weierstrass_p_inv(gr_ptr res, gr_srcptr z, gr_srcptr tau, gr_ctx_t ctx) { acb_elliptic_inv_p(res, z, tau, ACB_CTX_PREC(ctx)); return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE; }
int _gr_acb_weierstrass_zeta(gr_ptr res, gr_srcptr z, gr_srcptr tau, gr_ctx_t ctx) { acb_elliptic_zeta(res, z, tau, ACB_CTX_PREC(ctx)); return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE; }
int _gr_acb_weierstrass_sigma(gr_ptr res, gr_srcptr z, gr_srcptr tau, gr_ctx_t ctx) { acb_elliptic_sigma(res, z, tau, ACB_CTX_PREC(ctx)); return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE; }


int
_gr_acb_vec_dot(acb_t res, const acb_t initial, int subtract, acb_srcptr vec1, acb_srcptr vec2, slong len, gr_ctx_t ctx)
{
    acb_dot(res, initial, subtract, vec1, 1, vec2, 1, len, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_vec_dot_rev(acb_t res, const acb_t initial, int subtract, acb_srcptr vec1, acb_srcptr vec2, slong len, gr_ctx_t ctx)
{
    acb_dot(res, initial, subtract, vec1, 1, vec2 + len - 1, -1, len, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_poly_mullow(acb_ptr res,
    acb_srcptr poly1, slong len1,
    acb_srcptr poly2, slong len2, slong n, gr_ctx_t ctx)
{
    _acb_poly_mullow(res, poly1, len1, poly2, len2, n, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int _gr_acb_poly_taylor_shift(acb_ptr res, acb_srcptr poly, slong len, const acb_t c, gr_ctx_t ctx);

/* xxx */
static int
roots_accurate(acb_ptr roots, slong len, slong prec)
{
    slong i;

    for (i = 0; i < len; i++)
    {
        if (acb_rel_accuracy_bits(roots + i) < prec * 0.9)
            return 0;

        if (0)
        {
            if (mag_cmp_2exp_si(arb_radref(acb_realref(roots + i)), -prec * 0.9) > 0)
                return 0;

            if (mag_cmp_2exp_si(arb_radref(acb_imagref(roots + i)), -prec * 0.9) > 0)
                return 0;
        }
    }

    return 1;

}

/* hidden feature: also works with arb ctx */
int
_gr_acb_poly_roots(gr_vec_t roots, gr_vec_t mult, const gr_poly_t poly, int flags, gr_ctx_t ctx)
{
    slong prec, initial_prec, target_prec, isolated, maxiter, maxprec, deg, i;
    acb_ptr croots;
    acb_poly_t tmp;
    gr_ctx_t ZZ;
    int status = GR_UNABLE;
    int arb_roots;

    if (poly->length == 0)
        return GR_DOMAIN;

    deg = poly->length - 1;

    if (acb_contains_zero(((acb_srcptr) poly->coeffs) + deg))
        return GR_UNABLE;

    gr_ctx_init_fmpz(ZZ);
    croots = _acb_vec_init(deg);
    acb_poly_init(tmp);
    acb_poly_fit_length(tmp, deg + 1);
    _acb_poly_set_length(tmp, deg + 1);

    target_prec = ACB_CTX_PREC(ctx);
    initial_prec = 32;
    maxprec = 2 * target_prec + 64;

    if (ctx->which_ring == GR_CTX_RR_ARB)
        arb_roots = 1;
    else
        arb_roots = 0;

    for (prec = initial_prec; prec <= maxprec; prec *= 2)
    {
        maxiter = FLINT_MIN(2 * deg + 32, prec);
        for (i = 0; i <= deg; i++)
            _acb_vec_set_round(tmp->coeffs, poly->coeffs, deg + 1, prec);

        if (prec == initial_prec)
            isolated = acb_poly_find_roots(croots, tmp, NULL, maxiter, prec);
        else
            isolated = acb_poly_find_roots(croots, tmp, croots, maxiter, prec);

        if (isolated == deg)
        {
            if (arb_roots)
            {
                status = GR_UNABLE;

                if (roots_accurate(croots, deg, target_prec))
                {
                    if (acb_poly_validate_real_roots(croots, tmp, prec))
                        status = GR_SUCCESS;
                }
            }
            else
            {
                status = GR_SUCCESS;

                if (roots_accurate(croots, deg, target_prec))
                    break;
            }
        }
        else
        {
            status = GR_UNABLE;
        }
    }

    if (status == GR_SUCCESS)
    {
        _acb_vec_sort_pretty(croots, deg);

        if (arb_roots)
        {
            gr_vec_set_length(roots, 0, ctx);
            gr_vec_set_length(mult, 0, ZZ);

            for (i = 0; i < deg; i++)
            {
                if (arb_contains_zero(acb_imagref(croots + i)))
                {
                    fmpz one = 1;
                    arb_set_round(acb_realref(croots + i), acb_realref(croots + i), target_prec);
                    GR_MUST_SUCCEED(gr_vec_append(roots, acb_realref(croots + i), ctx));
                    GR_MUST_SUCCEED(gr_vec_append(mult, &one, ZZ));
                }
            }
        }
        else
        {
            gr_vec_set_length(roots, deg, ctx);
            gr_vec_set_length(mult, deg, ZZ);

            for (i = 0; i < deg; i++)
            {
                acb_set_round(((acb_ptr) roots->entries) + i, croots + i, target_prec);
                fmpz_one(((fmpz *) mult->entries) + i);
            }
        }
    }

    acb_poly_clear(tmp);
    _acb_vec_clear(croots, deg);
    gr_ctx_clear(ZZ);

    return status;
}

int
_gr_acb_poly_roots_other(gr_vec_t roots, gr_vec_t mult, const gr_poly_t poly, gr_ctx_t other_ctx, int flags, gr_ctx_t ctx)
{
    if (poly->length == 0)
        return GR_DOMAIN;

    if (other_ctx->which_ring == GR_CTX_CC_ACB)
        return _gr_acb_poly_roots(roots, mult, poly, flags, ctx);

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
                arb_fmpz_poly_complex_roots(croots, fac->p + i, 0, ACB_CTX_PREC(ctx));

                for (j = 0; j < deg2; j++)
                {
                    fmpz m2 = fac->exp[i];
                    GR_MUST_SUCCEED(gr_vec_append(roots, croots + j, ctx));
                    GR_MUST_SUCCEED(gr_vec_append(mult, &m2, ZZ));
                }

                _acb_vec_clear(croots, deg2);
            }

            fmpz_poly_factor_clear(fac);
        }

        gr_ctx_clear(ZZ);

        return status;
    }

    /* todo: first try conversion to fmpz? */
    {
        int status;
        gr_poly_t tmp;
        gr_poly_init(tmp, ctx);
        status = gr_poly_set_gr_poly_other(tmp, poly, other_ctx, ctx);
        if (status == GR_SUCCESS)
            status = _gr_acb_poly_roots(roots, mult, tmp, flags, ctx);
        gr_poly_clear(tmp, ctx);
        return status;
    }
}

int
_gr_acb_mat_mul(acb_mat_t res, const acb_mat_t x, const acb_mat_t y, gr_ctx_t ctx)
{
    acb_mat_mul(res, x, y, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_mat_det(acb_t res, const acb_mat_t x, gr_ctx_t ctx)
{
    acb_mat_det(res, x, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_mat_exp(acb_mat_t res, const acb_mat_t x, gr_ctx_t ctx)
{
    if (x->r != x->c)
        return GR_DOMAIN;

    acb_mat_exp(res, x, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_mat_diagonalization(gr_vec_t D, acb_mat_t L, acb_mat_t R, const acb_mat_t A, int flags, gr_ctx_t ctx)
{
    int status;
    slong n;
    acb_mat_t R_approx;

    n = acb_mat_nrows(A);

    if (n != acb_mat_ncols(A))
        return GR_DOMAIN;

    acb_mat_init(R_approx, n, n);

    gr_vec_set_length(D, n, ctx);

    acb_mat_approx_eig_qr(D->entries, NULL, R_approx, A, NULL, 0, ACB_CTX_PREC(ctx));
    status = acb_mat_eig_simple(D->entries, L, R, A, D->entries, R_approx, ACB_CTX_PREC(ctx)) ? GR_SUCCESS : GR_UNABLE;

    acb_mat_clear(R_approx);

    return status;
}


int _acb_methods_initialized = 0;

gr_static_method_table _acb_methods;

gr_method_tab_input _acb_methods_input[] =
{
    {GR_METHOD_CTX_WRITE,       (gr_funcptr) _gr_acb_ctx_write},
    {GR_METHOD_CTX_IS_RING,     (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_COMMUTATIVE_RING,
                                (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_INTEGRAL_DOMAIN,
                                (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_FIELD,    (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_UNIQUE_FACTORIZATION_DOMAIN,
                                (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_FINITE,
                                (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_FINITE_CHARACTERISTIC,
                                (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_ALGEBRAICALLY_CLOSED,
                                (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_ORDERED_RING,
                                (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_EXACT,    (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_CANONICAL,
                                (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_HAS_REAL_PREC, (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_SET_REAL_PREC, (gr_funcptr) _gr_acb_ctx_set_real_prec},
    {GR_METHOD_CTX_GET_REAL_PREC, (gr_funcptr) _gr_acb_ctx_get_real_prec},

    {GR_METHOD_INIT,            (gr_funcptr) _gr_acb_init},
    {GR_METHOD_CLEAR,           (gr_funcptr) _gr_acb_clear},
    {GR_METHOD_SWAP,            (gr_funcptr) _gr_acb_swap},
    {GR_METHOD_SET_SHALLOW,     (gr_funcptr) _gr_acb_set_shallow},
    {GR_METHOD_RANDTEST,        (gr_funcptr) _gr_acb_randtest},
    {GR_METHOD_WRITE,           (gr_funcptr) _gr_acb_write},
    {GR_METHOD_WRITE_N,         (gr_funcptr) _gr_acb_write_n},
    {GR_METHOD_ZERO,            (gr_funcptr) _gr_acb_zero},
    {GR_METHOD_ONE,             (gr_funcptr) _gr_acb_one},
    {GR_METHOD_IS_ZERO,         (gr_funcptr) _gr_acb_is_zero},
    {GR_METHOD_IS_ONE,          (gr_funcptr) _gr_acb_is_one},
    {GR_METHOD_IS_NEG_ONE,      (gr_funcptr) _gr_acb_is_neg_one},
    {GR_METHOD_EQUAL,           (gr_funcptr) _gr_acb_equal},
    {GR_METHOD_SET,             (gr_funcptr) _gr_acb_set},
    {GR_METHOD_SET_SI,          (gr_funcptr) _gr_acb_set_si},
    {GR_METHOD_SET_UI,          (gr_funcptr) _gr_acb_set_ui},
    {GR_METHOD_SET_FMPZ,        (gr_funcptr) _gr_acb_set_fmpz},
    {GR_METHOD_SET_FMPQ,        (gr_funcptr) _gr_acb_set_fmpq},
    {GR_METHOD_SET_OTHER,       (gr_funcptr) _gr_acb_set_other},
    {GR_METHOD_SET_D,           (gr_funcptr) _gr_acb_set_d},
    {GR_METHOD_GET_SI,          (gr_funcptr) _gr_acb_get_si},
    {GR_METHOD_GET_UI,          (gr_funcptr) _gr_acb_get_ui},
    {GR_METHOD_GET_FMPZ,        (gr_funcptr) _gr_acb_get_fmpz},
    {GR_METHOD_GET_D,           (gr_funcptr) _gr_acb_get_d},
    {GR_METHOD_NEG,             (gr_funcptr) _gr_acb_neg},
    {GR_METHOD_ADD,             (gr_funcptr) _gr_acb_add},
    {GR_METHOD_ADD_UI,          (gr_funcptr) _gr_acb_add_ui},
    {GR_METHOD_ADD_SI,          (gr_funcptr) _gr_acb_add_si},
    {GR_METHOD_ADD_FMPZ,        (gr_funcptr) _gr_acb_add_fmpz},
    {GR_METHOD_SUB,             (gr_funcptr) _gr_acb_sub},
    {GR_METHOD_SUB_UI,          (gr_funcptr) _gr_acb_sub_ui},
    {GR_METHOD_SUB_SI,          (gr_funcptr) _gr_acb_sub_si},
    {GR_METHOD_SUB_FMPZ,        (gr_funcptr) _gr_acb_sub_fmpz},
    {GR_METHOD_MUL,             (gr_funcptr) _gr_acb_mul},
    {GR_METHOD_MUL_UI,          (gr_funcptr) _gr_acb_mul_ui},
    {GR_METHOD_MUL_SI,          (gr_funcptr) _gr_acb_mul_si},
    {GR_METHOD_MUL_FMPZ,        (gr_funcptr) _gr_acb_mul_fmpz},
    {GR_METHOD_ADDMUL,          (gr_funcptr) _gr_acb_addmul},
    {GR_METHOD_SUBMUL,          (gr_funcptr) _gr_acb_submul},
    {GR_METHOD_MUL_TWO,         (gr_funcptr) _gr_acb_mul_two},
    {GR_METHOD_SQR,             (gr_funcptr) _gr_acb_sqr},
    {GR_METHOD_MUL_2EXP_SI,        (gr_funcptr) _gr_acb_mul_2exp_si},
    {GR_METHOD_MUL_2EXP_FMPZ,      (gr_funcptr) _gr_acb_mul_2exp_fmpz},
    {GR_METHOD_SET_FMPZ_2EXP_FMPZ, (gr_funcptr) _gr_acb_set_fmpz_2exp_fmpz},
    {GR_METHOD_GET_FMPZ_2EXP_FMPZ, (gr_funcptr) _gr_acb_get_fmpz_2exp_fmpz},
    {GR_METHOD_DIV,             (gr_funcptr) _gr_acb_div},
    {GR_METHOD_DIV_UI,          (gr_funcptr) _gr_acb_div_ui},
    {GR_METHOD_DIV_SI,          (gr_funcptr) _gr_acb_div_si},
    {GR_METHOD_DIV_FMPZ,        (gr_funcptr) _gr_acb_div_fmpz},
    {GR_METHOD_INV,             (gr_funcptr) _gr_acb_inv},
    {GR_METHOD_IS_INVERTIBLE,   (gr_funcptr) _gr_acb_is_invertible},
    {GR_METHOD_POW,             (gr_funcptr) _gr_acb_pow},
    {GR_METHOD_POW_UI,          (gr_funcptr) _gr_acb_pow_ui},
    {GR_METHOD_POW_SI,          (gr_funcptr) _gr_acb_pow_si},
    {GR_METHOD_POW_FMPZ,        (gr_funcptr) _gr_acb_pow_fmpz},
    {GR_METHOD_POW_FMPQ,        (gr_funcptr) _gr_acb_pow_fmpq},
    {GR_METHOD_IS_SQUARE,       (gr_funcptr) _gr_acb_is_square},
    {GR_METHOD_SQRT,            (gr_funcptr) _gr_acb_sqrt},
    {GR_METHOD_RSQRT,           (gr_funcptr) _gr_acb_rsqrt},
    {GR_METHOD_FLOOR,           (gr_funcptr) _gr_acb_floor},
    {GR_METHOD_CEIL,            (gr_funcptr) _gr_acb_ceil},
    {GR_METHOD_TRUNC,           (gr_funcptr) _gr_acb_trunc},
    {GR_METHOD_NINT,            (gr_funcptr) _gr_acb_nint},
    {GR_METHOD_I,               (gr_funcptr) _gr_acb_i},
    {GR_METHOD_ABS,             (gr_funcptr) _gr_acb_abs},
    {GR_METHOD_CONJ,            (gr_funcptr) _gr_acb_conj},
    {GR_METHOD_RE,              (gr_funcptr) _gr_acb_re},
    {GR_METHOD_IM,              (gr_funcptr) _gr_acb_im},
    {GR_METHOD_SGN,             (gr_funcptr) _gr_acb_sgn},
    {GR_METHOD_CSGN,            (gr_funcptr) _gr_acb_csgn},
    {GR_METHOD_ARG,             (gr_funcptr) _gr_acb_arg},
    {GR_METHOD_PI,              (gr_funcptr) _gr_acb_pi},
    {GR_METHOD_EXP,             (gr_funcptr) _gr_acb_exp},
    {GR_METHOD_EXPM1,           (gr_funcptr) _gr_acb_expm1},
    {GR_METHOD_EXP_PI_I,        (gr_funcptr) _gr_acb_exp_pi_i},
    {GR_METHOD_LOG,             (gr_funcptr) _gr_acb_log},
    {GR_METHOD_LOG1P,           (gr_funcptr) _gr_acb_log1p},
    {GR_METHOD_LOG_PI_I,        (gr_funcptr) _gr_acb_log_pi_i},
    {GR_METHOD_SIN,             (gr_funcptr) _gr_acb_sin},
    {GR_METHOD_COS,             (gr_funcptr) _gr_acb_cos},
    {GR_METHOD_SIN_COS,         (gr_funcptr) _gr_acb_sin_cos},
    {GR_METHOD_SIN_PI,          (gr_funcptr) _gr_acb_sin_pi},
    {GR_METHOD_COS_PI,          (gr_funcptr) _gr_acb_cos_pi},
    {GR_METHOD_SIN_COS_PI,      (gr_funcptr) _gr_acb_sin_cos_pi},
    {GR_METHOD_TAN,             (gr_funcptr) _gr_acb_tan},
    {GR_METHOD_COT,             (gr_funcptr) _gr_acb_cot},
    {GR_METHOD_SINC,            (gr_funcptr) _gr_acb_sinc},
    {GR_METHOD_SINC_PI,         (gr_funcptr) _gr_acb_sinc_pi},
    {GR_METHOD_SINH,            (gr_funcptr) _gr_acb_sinh},
    {GR_METHOD_COSH,            (gr_funcptr) _gr_acb_cosh},
    {GR_METHOD_SINH_COSH,       (gr_funcptr) _gr_acb_sinh_cosh},
    {GR_METHOD_TANH,            (gr_funcptr) _gr_acb_tanh},
    {GR_METHOD_COTH,            (gr_funcptr) _gr_acb_coth},
    {GR_METHOD_ASIN,            (gr_funcptr) _gr_acb_asin},
    {GR_METHOD_ACOS,            (gr_funcptr) _gr_acb_acos},
    {GR_METHOD_ATAN,            (gr_funcptr) _gr_acb_atan},
    {GR_METHOD_ASINH,           (gr_funcptr) _gr_acb_asinh},
    {GR_METHOD_ACOSH,           (gr_funcptr) _gr_acb_acosh},
    {GR_METHOD_ATANH,           (gr_funcptr) _gr_acb_atanh},
    {GR_METHOD_LAMBERTW,        (gr_funcptr) _gr_acb_lambertw},
    {GR_METHOD_LAMBERTW_FMPZ,   (gr_funcptr) _gr_acb_lambertw_fmpz},
    {GR_METHOD_FAC_UI,          (gr_funcptr) _gr_acb_fac_ui},
    {GR_METHOD_FAC_FMPZ,        (gr_funcptr) _gr_acb_fac_fmpz},
    {GR_METHOD_RISING_UI,       (gr_funcptr) _gr_acb_rising_ui},
    {GR_METHOD_RISING,          (gr_funcptr) _gr_acb_rising},
    {GR_METHOD_GAMMA,           (gr_funcptr) _gr_acb_gamma},
    {GR_METHOD_GAMMA_FMPZ,      (gr_funcptr) _gr_acb_gamma_fmpz},
    {GR_METHOD_GAMMA_FMPQ,      (gr_funcptr) _gr_acb_gamma_fmpq},
    {GR_METHOD_RGAMMA,          (gr_funcptr) _gr_acb_rgamma},
    {GR_METHOD_LGAMMA,          (gr_funcptr) _gr_acb_lgamma},
    {GR_METHOD_DIGAMMA,         (gr_funcptr) _gr_acb_digamma},
    {GR_METHOD_BARNES_G,        (gr_funcptr) _gr_acb_barnes_g},
    {GR_METHOD_LOG_BARNES_G,    (gr_funcptr) _gr_acb_log_barnes_g},
    {GR_METHOD_ERF,             (gr_funcptr) _gr_acb_erf},
    {GR_METHOD_ERFI,            (gr_funcptr) _gr_acb_erfi},
    {GR_METHOD_ERFC,            (gr_funcptr) _gr_acb_erfc},
    {GR_METHOD_ERFINV,          (gr_funcptr) _gr_acb_erfinv},
    {GR_METHOD_ERFCINV,         (gr_funcptr) _gr_acb_erfcinv},
    {GR_METHOD_FRESNEL_C,       (gr_funcptr) _gr_acb_fresnel_c},
    {GR_METHOD_FRESNEL_S,       (gr_funcptr) _gr_acb_fresnel_s},
    {GR_METHOD_FRESNEL,         (gr_funcptr) _gr_acb_fresnel},
    {GR_METHOD_GAMMA_UPPER,     (gr_funcptr) _gr_acb_gamma_upper},
    {GR_METHOD_GAMMA_LOWER,     (gr_funcptr) _gr_acb_gamma_lower},
    {GR_METHOD_BETA_LOWER,      (gr_funcptr) _gr_acb_beta_lower},
    {GR_METHOD_EXP_INTEGRAL,    (gr_funcptr) _gr_acb_exp_integral},
    {GR_METHOD_EXP_INTEGRAL_EI, (gr_funcptr) _gr_acb_exp_integral_ei},
    {GR_METHOD_SIN_INTEGRAL,    (gr_funcptr) _gr_acb_sin_integral},
    {GR_METHOD_COS_INTEGRAL,    (gr_funcptr) _gr_acb_cos_integral},
    {GR_METHOD_SINH_INTEGRAL,   (gr_funcptr) _gr_acb_sinh_integral},
    {GR_METHOD_COSH_INTEGRAL,   (gr_funcptr) _gr_acb_cosh_integral},
    {GR_METHOD_LOG_INTEGRAL,    (gr_funcptr) _gr_acb_log_integral},
    {GR_METHOD_DILOG,           (gr_funcptr) _gr_acb_dilog},
    {GR_METHOD_BESSEL_J,             (gr_funcptr) _gr_acb_bessel_j},
    {GR_METHOD_BESSEL_Y,             (gr_funcptr) _gr_acb_bessel_y},
    {GR_METHOD_BESSEL_I,             (gr_funcptr) _gr_acb_bessel_i},
    {GR_METHOD_BESSEL_K,             (gr_funcptr) _gr_acb_bessel_k},
    {GR_METHOD_BESSEL_J_Y,           (gr_funcptr) _gr_acb_bessel_j_y},
    {GR_METHOD_BESSEL_I_SCALED,      (gr_funcptr) _gr_acb_bessel_i_scaled},
    {GR_METHOD_BESSEL_K_SCALED,      (gr_funcptr) _gr_acb_bessel_k_scaled},
    {GR_METHOD_AIRY,                 (gr_funcptr) _gr_acb_airy},
    {GR_METHOD_AIRY_AI,              (gr_funcptr) _gr_acb_airy_ai},
    {GR_METHOD_AIRY_BI,              (gr_funcptr) _gr_acb_airy_bi},
    {GR_METHOD_AIRY_AI_PRIME,        (gr_funcptr) _gr_acb_airy_ai_prime},
    {GR_METHOD_AIRY_BI_PRIME,        (gr_funcptr) _gr_acb_airy_bi_prime},
    {GR_METHOD_AIRY_AI_ZERO,         (gr_funcptr) _gr_acb_airy_ai_zero},
    {GR_METHOD_AIRY_BI_ZERO,         (gr_funcptr) _gr_acb_airy_bi_zero},
    {GR_METHOD_AIRY_AI_PRIME_ZERO,   (gr_funcptr) _gr_acb_airy_ai_prime_zero},
    {GR_METHOD_AIRY_BI_PRIME_ZERO,   (gr_funcptr) _gr_acb_airy_bi_prime_zero},
    {GR_METHOD_COULOMB,              (gr_funcptr) _gr_acb_coulomb},
    {GR_METHOD_COULOMB_F,            (gr_funcptr) _gr_acb_coulomb_f},
    {GR_METHOD_COULOMB_G,            (gr_funcptr) _gr_acb_coulomb_g},
    {GR_METHOD_COULOMB_HNEG,         (gr_funcptr) _gr_acb_coulomb_hneg},
    {GR_METHOD_COULOMB_HPOS,         (gr_funcptr) _gr_acb_coulomb_hpos},
    {GR_METHOD_CHEBYSHEV_T,          (gr_funcptr) _gr_acb_chebyshev_t},
    {GR_METHOD_CHEBYSHEV_U,          (gr_funcptr) _gr_acb_chebyshev_u},
    {GR_METHOD_JACOBI_P,             (gr_funcptr) _gr_acb_jacobi_p},
    {GR_METHOD_GEGENBAUER_C,         (gr_funcptr) _gr_acb_gegenbauer_c},
    {GR_METHOD_LAGUERRE_L,           (gr_funcptr) _gr_acb_laguerre_l},
    {GR_METHOD_HERMITE_H,            (gr_funcptr) _gr_acb_hermite_h},
    {GR_METHOD_LEGENDRE_P,           (gr_funcptr) _gr_acb_legendre_p},
    {GR_METHOD_LEGENDRE_Q,           (gr_funcptr) _gr_acb_legendre_q},
    {GR_METHOD_LEGENDRE_P_ROOT_UI,   (gr_funcptr) _gr_acb_legendre_p_root_ui},
    {GR_METHOD_SPHERICAL_Y_SI,       (gr_funcptr) _gr_acb_spherical_y_si},
    {GR_METHOD_HYPGEOM_0F1,          (gr_funcptr) _gr_acb_hypgeom_0f1},
    {GR_METHOD_HYPGEOM_1F1,          (gr_funcptr) _gr_acb_hypgeom_1f1},
    {GR_METHOD_HYPGEOM_U,            (gr_funcptr) _gr_acb_hypgeom_u},
    {GR_METHOD_HYPGEOM_2F1,          (gr_funcptr) _gr_acb_hypgeom_2f1},
    {GR_METHOD_HYPGEOM_PFQ,          (gr_funcptr) _gr_acb_hypgeom_pfq},
    {GR_METHOD_ZETA,                 (gr_funcptr) _gr_acb_zeta},
    {GR_METHOD_HURWITZ_ZETA,         (gr_funcptr) _gr_acb_hurwitz_zeta},
    {GR_METHOD_POLYLOG,              (gr_funcptr) _gr_acb_polylog},
    {GR_METHOD_POLYGAMMA,            (gr_funcptr) _gr_acb_polygamma},
    {GR_METHOD_LERCH_PHI,            (gr_funcptr) _gr_acb_lerch_phi},
    {GR_METHOD_STIELTJES,            (gr_funcptr) _gr_acb_stieltjes},
    {GR_METHOD_DIRICHLET_ETA,        (gr_funcptr) _gr_acb_dirichlet_eta},
    {GR_METHOD_RIEMANN_XI,           (gr_funcptr) _gr_acb_riemann_xi},
    {GR_METHOD_ZETA_ZERO,            (gr_funcptr) _gr_acb_zeta_zero},
    {GR_METHOD_ZETA_ZERO_VEC,        (gr_funcptr) _gr_acb_zeta_zero_vec},
    {GR_METHOD_ZETA_NZEROS,          (gr_funcptr) _gr_acb_zeta_nzeros},

    {GR_METHOD_MODULAR_J,            (gr_funcptr) _gr_acb_modular_j},
    {GR_METHOD_MODULAR_DELTA,        (gr_funcptr) _gr_acb_modular_delta},
    {GR_METHOD_MODULAR_LAMBDA,       (gr_funcptr) _gr_acb_modular_lambda},
    {GR_METHOD_DEDEKIND_ETA,         (gr_funcptr) _gr_acb_dedekind_eta},
    {GR_METHOD_EISENSTEIN_G,         (gr_funcptr) _gr_acb_eisenstein_g},
    {GR_METHOD_EISENSTEIN_E,         (gr_funcptr) _gr_acb_eisenstein_e},
    {GR_METHOD_EISENSTEIN_G_VEC,     (gr_funcptr) _gr_acb_eisenstein_g_vec},

    {GR_METHOD_ELLIPTIC_K,           (gr_funcptr)  _gr_acb_elliptic_k},
    {GR_METHOD_ELLIPTIC_E,           (gr_funcptr)  _gr_acb_elliptic_e},
    {GR_METHOD_ELLIPTIC_PI,          (gr_funcptr)  _gr_acb_elliptic_pi},
    {GR_METHOD_ELLIPTIC_F,           (gr_funcptr)  _gr_acb_elliptic_f},
    {GR_METHOD_ELLIPTIC_E_INC,       (gr_funcptr)  _gr_acb_elliptic_e_inc},
    {GR_METHOD_ELLIPTIC_PI_INC,      (gr_funcptr)  _gr_acb_elliptic_pi_inc},

    {GR_METHOD_CARLSON_RC,           (gr_funcptr)  _gr_acb_carlson_rc},
    {GR_METHOD_CARLSON_RF,           (gr_funcptr)  _gr_acb_carlson_rf},
    {GR_METHOD_CARLSON_RG,           (gr_funcptr)  _gr_acb_carlson_rg},
    {GR_METHOD_CARLSON_RD,           (gr_funcptr)  _gr_acb_carlson_rd},
    {GR_METHOD_CARLSON_RJ,           (gr_funcptr)  _gr_acb_carlson_rj},

    {GR_METHOD_JACOBI_THETA,         (gr_funcptr)  _gr_acb_jacobi_theta},
    {GR_METHOD_JACOBI_THETA_1,       (gr_funcptr)  _gr_acb_jacobi_theta_1},
    {GR_METHOD_JACOBI_THETA_2,       (gr_funcptr)  _gr_acb_jacobi_theta_2},
    {GR_METHOD_JACOBI_THETA_3,       (gr_funcptr)  _gr_acb_jacobi_theta_3},
    {GR_METHOD_JACOBI_THETA_4,       (gr_funcptr)  _gr_acb_jacobi_theta_4},

    {GR_METHOD_ELLIPTIC_INVARIANTS,  (gr_funcptr)  _gr_acb_elliptic_invariants},
    {GR_METHOD_ELLIPTIC_ROOTS,       (gr_funcptr)  _gr_acb_elliptic_roots},

    {GR_METHOD_WEIERSTRASS_P,        (gr_funcptr)  _gr_acb_weierstrass_p},
    {GR_METHOD_WEIERSTRASS_P_PRIME,  (gr_funcptr)  _gr_acb_weierstrass_p_prime},
    {GR_METHOD_WEIERSTRASS_P_INV,    (gr_funcptr)  _gr_acb_weierstrass_p_inv},
    {GR_METHOD_WEIERSTRASS_ZETA,     (gr_funcptr)  _gr_acb_weierstrass_zeta},
    {GR_METHOD_WEIERSTRASS_SIGMA,    (gr_funcptr)  _gr_acb_weierstrass_sigma},

    {GR_METHOD_VEC_DOT,         (gr_funcptr) _gr_acb_vec_dot},
    {GR_METHOD_VEC_DOT_REV,     (gr_funcptr) _gr_acb_vec_dot_rev},

    {GR_METHOD_POLY_MULLOW,     (gr_funcptr) _gr_acb_poly_mullow},
    {GR_METHOD_POLY_TAYLOR_SHIFT,   (gr_funcptr) _gr_acb_poly_taylor_shift},
    {GR_METHOD_POLY_ROOTS,      (gr_funcptr) _gr_acb_poly_roots},
    {GR_METHOD_POLY_ROOTS_OTHER,(gr_funcptr) _gr_acb_poly_roots_other},
    {GR_METHOD_MAT_MUL,         (gr_funcptr) _gr_acb_mat_mul},
    {GR_METHOD_MAT_DET,         (gr_funcptr) _gr_acb_mat_det},
    {GR_METHOD_MAT_EXP,         (gr_funcptr) _gr_acb_mat_exp},
    {GR_METHOD_MAT_DIAGONALIZATION,     (gr_funcptr) _gr_acb_mat_diagonalization},
    {0,                         (gr_funcptr) NULL},
};

void
gr_ctx_init_complex_acb(gr_ctx_t ctx, slong prec)
{
    ctx->which_ring = GR_CTX_CC_ACB;
    ctx->sizeof_elem = sizeof(acb_struct);
    ctx->size_limit = WORD_MAX;

    ACB_CTX_PREC(ctx) = FLINT_MAX(2, FLINT_MIN(prec, WORD_MAX / 8));
    ACB_CTX_PREC(ctx) = prec;

    ctx->methods = _acb_methods;

    if (!_acb_methods_initialized)
    {
        gr_method_tab_init(_acb_methods, _acb_methods_input);
        _acb_methods_initialized = 1;
    }
}
