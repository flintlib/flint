/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpzi.h"
#include "ca.h"
#include "ca_mat.h"
#include "ca_poly.h"
#include "fexpr.h"
#include "gr.h"
#include "gr_generic.h"
#include "gr_vec.h"
#include "gr_poly.h"
#include "gr_special.h"

#define GR_CA_CTX(ring_ctx) ((ca_ctx_struct *)(GR_CTX_DATA_AS_PTR(ring_ctx)))

int
_gr_ca_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    if (ctx->which_ring == GR_CTX_RR_CA)
        gr_stream_write(out, "Real numbers (ca)");
    else if (ctx->which_ring == GR_CTX_CC_CA)
        gr_stream_write(out, "Complex numbers (ca)");
    else if (ctx->which_ring == GR_CTX_REAL_ALGEBRAIC_CA)
        gr_stream_write(out, "Real algebraic numbers (ca)");
    else if (ctx->which_ring == GR_CTX_COMPLEX_ALGEBRAIC_CA)
        gr_stream_write(out, "Complex algebraic numbers (ca)");
    else if (ctx->which_ring == GR_CTX_COMPLEX_EXTENDED_CA)
        gr_stream_write(out, "Complex numbers + extended values (ca)");
    return GR_SUCCESS;
}

void
gr_ctx_ca_set_option(gr_ctx_t ctx, slong option, slong value)
{
    ca_ctx_set_option(GR_CA_CTX(ctx), option, value);
}

slong
gr_ctx_ca_get_option(gr_ctx_t ctx, slong option)
{
    return ca_ctx_get_option(GR_CA_CTX(ctx), option);
}

void
_gr_ca_init(ca_t x, gr_ctx_t ctx)
{
    ca_init(x, GR_CA_CTX(ctx));
}

void
_gr_ca_clear(ca_t x, gr_ctx_t ctx)
{
    ca_clear(x, GR_CA_CTX(ctx));
}

void
_gr_ca_swap(ca_t x, ca_t y, gr_ctx_t ctx)
{
    ca_t t;
    *t = *x;
    *x = *y;
    *y = *t;
}

void
_gr_ca_set_shallow(ca_t res, const ca_t x, const gr_ctx_t ctx)
{
    *res = *x;
}

/* todo: limits */
/* todo: faster real/algebraic constructions */
int
_gr_ca_randtest(ca_t res, flint_rand_t state, gr_ctx_t ctx)
{
    if (ctx->which_ring == GR_CTX_COMPLEX_EXTENDED_CA)
    {
        ca_randtest_special(res, state, 2, 10, GR_CA_CTX(ctx));
    }
    else
    {
        ca_randtest(res, state, 2, 10, GR_CA_CTX(ctx));

        if (ctx->which_ring == GR_CTX_RR_CA)
        {
            if (ca_check_is_real(res, GR_CA_CTX(ctx)) != T_TRUE)
            {
                ca_randtest_rational(res, state, 10, GR_CA_CTX(ctx));
            }
        }
        else if (ctx->which_ring == GR_CTX_REAL_ALGEBRAIC_CA)
        {
            if (ca_check_is_real(res, GR_CA_CTX(ctx)) != T_TRUE ||
                ca_check_is_algebraic(res, GR_CA_CTX(ctx)) != T_TRUE)
            {
                ca_randtest_rational(res, state, 10, GR_CA_CTX(ctx));
            }
        }
        else if (ctx->which_ring == GR_CTX_COMPLEX_ALGEBRAIC_CA)
        {
            if (ca_check_is_algebraic(res, GR_CA_CTX(ctx)) != T_TRUE)
            {
                ca_randtest_rational(res, state, 10, GR_CA_CTX(ctx));
            }
        }
    }

    return GR_SUCCESS;
}

/* todo */
int
_gr_ca_write(gr_stream_t out, const ca_t x, gr_ctx_t ctx)
{
    gr_stream_write_free(out, ca_get_str(x, GR_CA_CTX(ctx)));
    return GR_SUCCESS;
}

int
_gr_ca_zero(ca_t x, gr_ctx_t ctx)
{
    ca_zero(x, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_one(ca_t x, gr_ctx_t ctx)
{
    ca_one(x, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_set_si(ca_t res, slong v, gr_ctx_t ctx)
{
    ca_set_si(res, v, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_set_ui(ca_t res, ulong v, gr_ctx_t ctx)
{
    ca_set_ui(res, v, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_set_fmpz(ca_t res, const fmpz_t v, gr_ctx_t ctx)
{
    ca_set_fmpz(res, v, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_set_fmpq(ca_t res, const fmpq_t v, gr_ctx_t ctx)
{
    ca_set_fmpq(res, v, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

void
ca_set_fmpzi(ca_t res, const fmpzi_t x, ca_ctx_t ctx)
{
    if (fmpz_is_zero(fmpzi_imagref(x)))
    {
        ca_set_fmpz(res, fmpzi_realref(x), ctx);
    }
    else
    {
        ca_i(res, ctx);
        ca_mul_fmpz(res, res, fmpzi_imagref(x), ctx);
        ca_add_fmpz(res, res, fmpzi_realref(x), ctx);
    }
}

int
_gr_ca_set_other(ca_t res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
{
    slong target = ctx->which_ring;

    switch (x_ctx->which_ring)
    {
        case GR_CTX_FMPZ:
            ca_set_fmpz(res, x, GR_CA_CTX(ctx));
            return GR_SUCCESS;

        case GR_CTX_FMPQ:
            ca_set_fmpq(res, x, GR_CA_CTX(ctx));
            return GR_SUCCESS;

        case GR_CTX_FMPZI:
            if (target == GR_CTX_CC_CA || target == GR_CTX_COMPLEX_ALGEBRAIC_CA || target == GR_CTX_COMPLEX_EXTENDED_CA ||
                fmpz_is_zero(fmpzi_imagref((const fmpzi_struct *) x)))
            {
                ca_set_fmpzi(res, x, GR_CA_CTX(ctx));
                return GR_SUCCESS;
            }
            else
            {
                return GR_DOMAIN;
            }

        case GR_CTX_REAL_ALGEBRAIC_QQBAR:
            ca_set_qqbar(res, x, GR_CA_CTX(ctx));
            return GR_SUCCESS;

        case GR_CTX_COMPLEX_ALGEBRAIC_QQBAR:
            if (target == GR_CTX_CC_CA ||
                target == GR_CTX_COMPLEX_ALGEBRAIC_CA ||
                target == GR_CTX_COMPLEX_EXTENDED_CA ||
                qqbar_is_real(x))
            {
                ca_set_qqbar(res, x, GR_CA_CTX(ctx));
                return GR_SUCCESS;
            }
            else
            {
                return GR_DOMAIN;
            }

        case GR_CTX_REAL_ALGEBRAIC_CA:
            ca_transfer(res, GR_CA_CTX(ctx), x, GR_CA_CTX(x_ctx));
            return GR_SUCCESS;

        case GR_CTX_COMPLEX_ALGEBRAIC_CA:
            {
                truth_t ok = T_UNKNOWN;

                if (target == GR_CTX_CC_CA || target == GR_CTX_COMPLEX_ALGEBRAIC_CA || target == GR_CTX_COMPLEX_EXTENDED_CA)
                {
                    ok = T_TRUE;
                }
                else if (target == GR_CTX_RR_CA)
                {
                    ok = ca_check_is_real(x, GR_CA_CTX(x_ctx));
                }
                else if (target == GR_CTX_REAL_ALGEBRAIC_CA)
                {
                    ok = truth_and(ca_check_is_algebraic(x, GR_CA_CTX(x_ctx)), ca_check_is_real(x, GR_CA_CTX(x_ctx)));
                }

                if (ok == T_TRUE)
                {
                    ca_transfer(res, GR_CA_CTX(ctx), x, GR_CA_CTX(x_ctx));
                    return GR_SUCCESS;
                }

                return (ok == T_FALSE) ? GR_DOMAIN : GR_UNABLE;
            }

        case GR_CTX_RR_CA:
            {
                truth_t ok = T_UNKNOWN;

                if (target == GR_CTX_RR_CA || target == GR_CTX_CC_CA || target == GR_CTX_COMPLEX_EXTENDED_CA)
                {
                    ok = T_TRUE;
                }
                else if (target == GR_CTX_REAL_ALGEBRAIC_CA || target == GR_CTX_COMPLEX_ALGEBRAIC_CA)
                {
                    ok = ca_check_is_algebraic(x, GR_CA_CTX(x_ctx));
                }

                if (ok == T_TRUE)
                {
                    ca_transfer(res, GR_CA_CTX(ctx), x, GR_CA_CTX(x_ctx));
                    return GR_SUCCESS;
                }

                return (ok == T_FALSE) ? GR_DOMAIN : GR_UNABLE;
            }

        case GR_CTX_CC_CA:
            {
                truth_t ok = T_UNKNOWN;

                if (target == GR_CTX_CC_CA || target == GR_CTX_COMPLEX_EXTENDED_CA)
                {
                    ok = T_TRUE;
                }
                else if (target == GR_CTX_RR_CA)
                {
                    ok = ca_check_is_real(x, GR_CA_CTX(x_ctx));
                }
                else if (target == GR_CTX_COMPLEX_ALGEBRAIC_CA)
                {
                    ok = ca_check_is_algebraic(x, GR_CA_CTX(x_ctx));
                }
                else if (target == GR_CTX_REAL_ALGEBRAIC_CA)
                {
                    ok = truth_and(ca_check_is_algebraic(x, GR_CA_CTX(x_ctx)), ca_check_is_real(x, GR_CA_CTX(x_ctx)));
                }

                if (ok == T_TRUE)
                {
                    ca_transfer(res, GR_CA_CTX(ctx), x, GR_CA_CTX(x_ctx));
                    return GR_SUCCESS;
                }

                return (ok == T_FALSE) ? GR_DOMAIN : GR_UNABLE;
            }

        case GR_CTX_COMPLEX_EXTENDED_CA:
            {
                truth_t ok = T_UNKNOWN;

                if (target == GR_CTX_COMPLEX_EXTENDED_CA)
                {
                    ok = T_TRUE;
                }
                else
                {
                    if (ca_check_is_undefined(x, GR_CA_CTX(x_ctx)) == T_TRUE ||
                         ca_check_is_infinity(x, GR_CA_CTX(x_ctx)) == T_TRUE)
                    {
                        ok = T_FALSE;
                    }
                    else if (ca_is_unknown(x, GR_CA_CTX(x_ctx)))
                    {
                        ok = T_UNKNOWN;
                    }
                    else
                    {
                        if (target == GR_CTX_RR_CA)
                        {
                            ok = ca_check_is_real(x, GR_CA_CTX(x_ctx));
                        }
                        else if (target == GR_CTX_COMPLEX_ALGEBRAIC_CA)
                        {
                            ok = ca_check_is_algebraic(x, GR_CA_CTX(x_ctx));
                        }
                        else if (target == GR_CTX_REAL_ALGEBRAIC_CA)
                        {
                            ok = truth_and(ca_check_is_algebraic(x, GR_CA_CTX(x_ctx)), ca_check_is_real(x, GR_CA_CTX(x_ctx)));
                        }
                        else
                        {
                            ok = T_TRUE;
                        }
                    }
                }

                if (ok == T_TRUE)
                {
                    ca_transfer(res, GR_CA_CTX(ctx), x, GR_CA_CTX(x_ctx));
                    return GR_SUCCESS;
                }

                return (ok == T_FALSE) ? GR_DOMAIN : GR_UNABLE;
            }

    }

    return gr_generic_set_other(res, x, x_ctx, ctx);
}

int
_gr_ca_get_arb_with_prec(arb_t res, gr_srcptr x, gr_ctx_t x_ctx, slong prec)
{
    int status = GR_UNABLE;
    truth_t ok;

    acb_t t;
    acb_init(t);

    /* todo: when to use accurate_parts? */
    ca_get_acb(t, x, prec, GR_CA_CTX(x_ctx));

    if (x_ctx->which_ring == GR_CTX_RR_CA || x_ctx->which_ring == GR_CTX_REAL_ALGEBRAIC_CA ||
        (arb_is_zero(acb_imagref(t)) && arb_is_finite(acb_realref(t))))
    {
        status = GR_SUCCESS;
    }
    else
    {
        ok = ca_check_is_real(x, GR_CA_CTX(x_ctx));

        if (ok == T_TRUE)
            status = GR_SUCCESS;
        else
            status = (ok == T_FALSE) ? GR_DOMAIN : GR_UNABLE;
    }

    if (status == GR_SUCCESS)
        arb_set_round(res, acb_realref(t), prec);

    acb_clear(t);
    return status;
}

int
_gr_ca_get_acb_with_prec(acb_t res, gr_srcptr x, gr_ctx_t x_ctx, slong prec)
{
    if (x_ctx->which_ring == GR_CTX_COMPLEX_EXTENDED_CA)
    {
        if (ca_check_is_undefined(x, GR_CA_CTX(x_ctx)) == T_TRUE ||
             ca_check_is_infinity(x, GR_CA_CTX(x_ctx)) == T_TRUE)
            return GR_DOMAIN;

        if (ca_is_unknown(x, GR_CA_CTX(x_ctx)))
            return GR_UNABLE;
    }

    /* todo: when to use accurate_parts? */
    ca_get_acb(res, x, prec, GR_CA_CTX(x_ctx));
    acb_set_round(res, res, prec);
    return GR_SUCCESS;
}

int
_gr_ca_get_fmpz(fmpz_t res, const ca_t x, gr_ctx_t ctx)
{
    truth_t integer;

    /* todo: avoid duplicate computation */
    integer = ca_check_is_integer(x, GR_CA_CTX(ctx));

    if (integer == T_TRUE)
        return ca_get_fmpz(res, x, GR_CA_CTX(ctx)) ? GR_SUCCESS : GR_UNABLE;
    else if (integer == T_FALSE)
        return GR_DOMAIN;
    else
        return GR_UNABLE;
}

int
_gr_ca_get_si(slong * res, const ca_t x, gr_ctx_t ctx)
{
    fmpz_t n;
    int status;

    fmpz_init(n);
    status = _gr_ca_get_fmpz(n, x, ctx);

    if (status == GR_SUCCESS)
    {
        if (fmpz_fits_si(n))
            *res = fmpz_get_si(n);
        else
            status = GR_DOMAIN;
    }

    fmpz_clear(n);
    return status;
}

int
_gr_ca_get_ui(ulong * res, const ca_t x, gr_ctx_t ctx)
{
    fmpz_t n;
    int status;

    fmpz_init(n);
    status = _gr_ca_get_fmpz(n, x, ctx);

    if (status == GR_SUCCESS)
    {
        if (fmpz_sgn(n) >= 0 && fmpz_cmp_ui(n, UWORD_MAX) <= 0)
            *res = fmpz_get_ui(n);
        else
            status = GR_DOMAIN;
    }

    fmpz_clear(n);
    return status;
}

int
_gr_ca_get_d(double * res, gr_srcptr x, gr_ctx_t ctx)
{
    arb_t t;
    int status = GR_UNABLE;

    arb_init(t);
    status = _gr_ca_get_arb_with_prec(t, x, ctx, 64);
    if (status == GR_SUCCESS)
        *res = arf_get_d(arb_midref(t), ARF_RND_NEAR);
    arb_clear(t);
    return status;
}

int
_gr_ca_get_fmpq(fmpq_t res, const ca_t x, gr_ctx_t ctx)
{
    truth_t rational;

    /* todo: avoid duplicate computation */
    rational = ca_check_is_rational(x, GR_CA_CTX(ctx));

    if (rational == T_TRUE)
        return ca_get_fmpq(res, x, GR_CA_CTX(ctx)) ? GR_SUCCESS : GR_UNABLE;
    else if (rational == T_FALSE)
        return GR_DOMAIN;
    else
        return GR_UNABLE;
}

int
_gr_ca_get_fexpr(fexpr_t res, const ca_t x, gr_ctx_t ctx)
{
    ca_get_fexpr(res, x, 0, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_get_fexpr_serialize(fexpr_t res, const ca_t x, gr_ctx_t ctx)
{
    ca_get_fexpr(res, x, CA_FEXPR_SERIALIZATION, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}


truth_t
_gr_ca_is_zero(const ca_t x, gr_ctx_t ctx)
{
    return ca_check_is_zero(x, GR_CA_CTX(ctx));
}

truth_t
_gr_ca_is_one(const ca_t x, gr_ctx_t ctx)
{
    return ca_check_is_one(x, GR_CA_CTX(ctx));
}

truth_t
_gr_ca_is_neg_one(const ca_t x, gr_ctx_t ctx)
{
    return ca_check_is_neg_one(x, GR_CA_CTX(ctx));
}

truth_t
_gr_ca_equal(const ca_t x, const ca_t y, gr_ctx_t ctx)
{
    return ca_check_equal(x, y, GR_CA_CTX(ctx));
}

int
_gr_ca_set(ca_t res, const ca_t x, gr_ctx_t ctx)
{
    ca_set(res, x, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_neg(ca_t res, const ca_t x, gr_ctx_t ctx)
{
    ca_neg(res, x, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_add(ca_t res, const ca_t x, const ca_t y, gr_ctx_t ctx)
{
    ca_add(res, x, y, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_add_si(ca_t res, const ca_t x, slong y, gr_ctx_t ctx)
{
    ca_add_si(res, x, y, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_add_ui(ca_t res, const ca_t x, ulong y, gr_ctx_t ctx)
{
    ca_add_ui(res, x, y, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_add_fmpz(ca_t res, const ca_t x, const fmpz_t y, gr_ctx_t ctx)
{
    ca_add_fmpz(res, x, y, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_add_fmpq(ca_t res, const ca_t x, const fmpq_t y, gr_ctx_t ctx)
{
    ca_add_fmpq(res, x, y, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_sub(ca_t res, const ca_t x, const ca_t y, gr_ctx_t ctx)
{
    ca_sub(res, x, y, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_sub_si(ca_t res, const ca_t x, slong y, gr_ctx_t ctx)
{
    ca_sub_si(res, x, y, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_sub_ui(ca_t res, const ca_t x, ulong y, gr_ctx_t ctx)
{
    ca_sub_ui(res, x, y, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_sub_fmpz(ca_t res, const ca_t x, const fmpz_t y, gr_ctx_t ctx)
{
    ca_sub_fmpz(res, x, y, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_sub_fmpq(ca_t res, const ca_t x, const fmpq_t y, gr_ctx_t ctx)
{
    ca_sub_fmpq(res, x, y, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_mul(ca_t res, const ca_t x, const ca_t y, gr_ctx_t ctx)
{
    ca_mul(res, x, y, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_mul_si(ca_t res, const ca_t x, slong y, gr_ctx_t ctx)
{
    ca_mul_si(res, x, y, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_mul_ui(ca_t res, const ca_t x, ulong y, gr_ctx_t ctx)
{
    ca_mul_ui(res, x, y, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_mul_fmpz(ca_t res, const ca_t x, const fmpz_t y, gr_ctx_t ctx)
{
    ca_mul_fmpz(res, x, y, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_mul_fmpq(ca_t res, const ca_t x, const fmpq_t y, gr_ctx_t ctx)
{
    ca_mul_fmpq(res, x, y, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

static int
handle_possible_special_value(ca_t res, gr_ctx_t ctx)
{
    if (ctx->which_ring == GR_CTX_COMPLEX_EXTENDED_CA)
        return GR_SUCCESS;

    if (ca_is_unknown(res, GR_CA_CTX(ctx)))
        return GR_UNABLE;

    if (ca_is_special(res, GR_CA_CTX(ctx)))
    {
        ca_unknown(res, GR_CA_CTX(ctx));   /* don't output an invalid element */
        return GR_DOMAIN;
    }

    return GR_SUCCESS;
}

#define DEF_SPECIAL(name) \
int \
_gr_ca_ ## name(ca_t res, gr_ctx_t ctx) \
{ \
    if (ctx->which_ring != GR_CTX_COMPLEX_EXTENDED_CA) \
        return GR_DOMAIN; \
    ca_ ## name(res, GR_CA_CTX(ctx)); \
    return GR_SUCCESS; \
}

DEF_SPECIAL(pos_inf)
DEF_SPECIAL(neg_inf)
DEF_SPECIAL(uinf)
DEF_SPECIAL(undefined)
DEF_SPECIAL(unknown)

int
_gr_ca_inv(ca_t res, const ca_t x, gr_ctx_t ctx)
{
    ca_inv(res, x, GR_CA_CTX(ctx));
    return handle_possible_special_value(res, ctx);
}

int
_gr_ca_div(ca_t res, const ca_t x, const ca_t y, gr_ctx_t ctx)
{
    ca_div(res, x, y, GR_CA_CTX(ctx));
    return handle_possible_special_value(res, ctx);
}

int
_gr_ca_div_si(ca_t res, const ca_t x, slong y, gr_ctx_t ctx)
{
    ca_div_si(res, x, y, GR_CA_CTX(ctx));
    return handle_possible_special_value(res, ctx);
}

int
_gr_ca_div_ui(ca_t res, const ca_t x, ulong y, gr_ctx_t ctx)
{
    ca_div_ui(res, x, y, GR_CA_CTX(ctx));
    return handle_possible_special_value(res, ctx);
}

int
_gr_ca_div_fmpz(ca_t res, const ca_t x, const fmpz_t y, gr_ctx_t ctx)
{
    ca_div_fmpz(res, x, y, GR_CA_CTX(ctx));
    return handle_possible_special_value(res, ctx);
}

int
_gr_ca_div_fmpq(ca_t res, const ca_t x, const fmpq_t y, gr_ctx_t ctx)
{
    ca_div_fmpq(res, x, y, GR_CA_CTX(ctx));
    return handle_possible_special_value(res, ctx);
}

truth_t
_gr_ca_is_invertible(const ca_t x, gr_ctx_t ctx)
{
    return truth_not(ca_check_is_zero(x, GR_CA_CTX(ctx)));
}

int
_gr_ca_pow_ui(ca_t res, const ca_t x, ulong exp, gr_ctx_t ctx)
{
    ca_pow_ui(res, x, exp, GR_CA_CTX(ctx));
    return handle_possible_special_value(res, ctx);
}

int
_gr_ca_pow_si(ca_t res, const ca_t x, slong exp, gr_ctx_t ctx)
{
    ca_pow_si(res, x, exp, GR_CA_CTX(ctx));
    return handle_possible_special_value(res, ctx);
}

int
_gr_ca_pow_fmpz(ca_t res, const ca_t x, const fmpz_t exp, gr_ctx_t ctx)
{
    ca_pow_fmpz(res, x, exp, GR_CA_CTX(ctx));
    return handle_possible_special_value(res, ctx);
}

int
_gr_ca_pow_fmpq(ca_t res, const ca_t x, const fmpq_t exp, gr_ctx_t ctx)
{
    ca_pow_fmpq(res, x, exp, GR_CA_CTX(ctx));

    if (ctx->which_ring == GR_CTX_RR_CA || ctx->which_ring == GR_CTX_REAL_ALGEBRAIC_CA)
    {
        truth_t real;

        real = ca_check_is_real(res, GR_CA_CTX(ctx));

        if (real == T_UNKNOWN)
            return GR_UNABLE;

        if (real == T_FALSE)
            return GR_DOMAIN;
    }

    return handle_possible_special_value(res, ctx);
}

int
_gr_ca_pow(ca_t res, const ca_t x, const ca_t exp, gr_ctx_t ctx)
{
    ca_pow(res, x, exp, GR_CA_CTX(ctx));

    if (ctx->which_ring == GR_CTX_RR_CA || ctx->which_ring == GR_CTX_REAL_ALGEBRAIC_CA)
    {
        truth_t real;

        real = ca_check_is_real(res, GR_CA_CTX(ctx));

        if (real == T_UNKNOWN)
            return GR_UNABLE;

        if (real == T_FALSE)
            return GR_DOMAIN;
    }

    if (ctx->which_ring == GR_CTX_REAL_ALGEBRAIC_CA || ctx->which_ring == GR_CTX_COMPLEX_ALGEBRAIC_CA)
    {
        truth_t algebraic;

        algebraic = ca_check_is_algebraic(res, GR_CA_CTX(ctx));

        if (algebraic == T_UNKNOWN)
            return GR_UNABLE;

        if (algebraic == T_FALSE)
            return GR_DOMAIN;
    }

    return handle_possible_special_value(res, ctx);
}

truth_t
_gr_ca_is_square(const ca_t x, gr_ctx_t ctx)
{
    if (ctx->which_ring == GR_CTX_RR_CA || ctx->which_ring == GR_CTX_REAL_ALGEBRAIC_CA)
    {
        return truth_not(ca_check_is_negative_real(x, GR_CA_CTX(ctx)));
    }
    else
    {
        return T_TRUE;
    }
}

int
_gr_ca_sqrt(ca_t res, const ca_t x, gr_ctx_t ctx)
{
    ca_sqrt(res, x, GR_CA_CTX(ctx));

    if (ctx->which_ring == GR_CTX_RR_CA || ctx->which_ring == GR_CTX_REAL_ALGEBRAIC_CA)
    {
        truth_t real;

        real = ca_check_is_real(res, GR_CA_CTX(ctx));

        if (real == T_UNKNOWN)
            return GR_UNABLE;

        if (real == T_FALSE)
            return GR_DOMAIN;
    }

    return handle_possible_special_value(res, ctx);
}

int
_gr_ca_rsqrt(ca_t res, const ca_t x, gr_ctx_t ctx)
{
    ca_sqrt(res, x, GR_CA_CTX(ctx));
    ca_inv(res, res, GR_CA_CTX(ctx));

    if (ctx->which_ring == GR_CTX_RR_CA || ctx->which_ring == GR_CTX_REAL_ALGEBRAIC_CA)
    {
        truth_t real;

        real = ca_check_is_real(res, GR_CA_CTX(ctx));

        if (real == T_UNKNOWN)
            return GR_UNABLE;

        if (real == T_FALSE)
            return GR_DOMAIN;
    }

    return handle_possible_special_value(res, ctx);
}

int
_gr_ca_floor(ca_t res, const ca_t x, gr_ctx_t ctx)
{
    ca_floor(res, x, GR_CA_CTX(ctx));

    if (ca_is_unknown(res, GR_CA_CTX(ctx)))
        return GR_UNABLE;

    return GR_SUCCESS;
}

int
_gr_ca_ceil(ca_t res, const ca_t x, gr_ctx_t ctx)
{
    ca_ceil(res, x, GR_CA_CTX(ctx));

    if (ca_is_unknown(res, GR_CA_CTX(ctx)))
        return GR_UNABLE;

    return GR_SUCCESS;
}

/* todo: trunc, nint in calcium */

int
_gr_ca_trunc(ca_t res, const ca_t x, gr_ctx_t ctx)
{
    acb_t t;
    int status;

    acb_init(t);

    ca_get_acb(t, x, 64, GR_CA_CTX(ctx));

    if (arf_cmpabs_2exp_si(arb_midref(acb_realref(t)), -1) < 0 && mag_cmp_2exp_si(arb_radref(acb_realref(t)), -1) < 0)
    {
        ca_zero(res, GR_CA_CTX(ctx));
        status = GR_SUCCESS;
    }
    else if (arb_is_positive(acb_realref(t)))
    {
        status = _gr_ca_floor(res, x, ctx);
    }
    else if (arb_is_negative(acb_realref(t)))
    {
        status = _gr_ca_ceil(res, x, ctx);
    }
    else
    {
        status = GR_UNABLE;
    }

    acb_clear(t);

    return status;
}

/* todo: fast numerical path */
int
_gr_ca_nint(ca_t res, const ca_t x, gr_ctx_t ctx)
{
    if (ca_check_is_integer(x, GR_CA_CTX(ctx)) == T_TRUE)
    {
        ca_set(res, x, GR_CA_CTX(ctx));
        return GR_SUCCESS;
    }
    else
    {
        ca_t t;
        truth_t integer;
        int status = GR_SUCCESS;

        ca_init(t, GR_CA_CTX(ctx));

        ca_set_d(t, 0.5, GR_CA_CTX(ctx));
        ca_add(t, x, t, GR_CA_CTX(ctx));
        ca_re(t, t, GR_CA_CTX(ctx));
        ca_floor(res, t, GR_CA_CTX(ctx));

        integer = ca_check_is_integer(t, GR_CA_CTX(ctx));

        if (integer == T_TRUE)
        {
            fmpz_t m;
            fmpz_init(m);
            if (ca_get_fmpz(m, t, GR_CA_CTX(ctx)))
            {
                if (fmpz_is_odd(m))
                    ca_sub_ui(res, res, 1, GR_CA_CTX(ctx));
            }
            else
            {
                status = GR_UNABLE;
            }
            fmpz_clear(m);
        }
        else if (integer == T_UNKNOWN)
        {
            status = GR_UNABLE;
        }

        ca_clear(t, GR_CA_CTX(ctx));
        return status;
    }
}

int
_gr_ca_i(ca_t res, gr_ctx_t ctx)
{
    if (ctx->which_ring == GR_CTX_REAL_ALGEBRAIC_CA ||
        ctx->which_ring == GR_CTX_RR_CA)
        return GR_DOMAIN;

    ca_i(res, GR_CA_CTX(ctx));

    if (ca_is_unknown(res, GR_CA_CTX(ctx)))
        return GR_UNABLE;

    return GR_SUCCESS;
}

int
_gr_ca_abs(ca_t res, const ca_t x, gr_ctx_t ctx)
{
    ca_abs(res, x, GR_CA_CTX(ctx));

    if (ca_is_unknown(res, GR_CA_CTX(ctx)))
        return GR_UNABLE;

    return GR_SUCCESS;
}

int
_gr_ca_conj(ca_t res, const ca_t x, gr_ctx_t ctx)
{
    ca_conj(res, x, GR_CA_CTX(ctx));

    if (ca_is_unknown(res, GR_CA_CTX(ctx)))
        return GR_UNABLE;

    return GR_SUCCESS;
}

/* todo: exploit when we know that the field is real */
int
_gr_ca_re(ca_t res, const ca_t x, gr_ctx_t ctx)
{
    ca_re(res, x, GR_CA_CTX(ctx));

    if (ca_is_unknown(res, GR_CA_CTX(ctx)))
        return GR_UNABLE;

    return GR_SUCCESS;
}

/* todo: exploit when we know that the field is real */
int
_gr_ca_im(ca_t res, const ca_t x, gr_ctx_t ctx)
{
    ca_im(res, x, GR_CA_CTX(ctx));

    if (ca_is_unknown(res, GR_CA_CTX(ctx)))
        return GR_UNABLE;

    return GR_SUCCESS;
}

int
_gr_ca_sgn(ca_t res, const ca_t x, gr_ctx_t ctx)
{
    ca_sgn(res, x, GR_CA_CTX(ctx));

    if (ca_is_unknown(res, GR_CA_CTX(ctx)))
        return GR_UNABLE;

    return GR_SUCCESS;
}

int
_gr_ca_csgn(ca_t res, const ca_t x, gr_ctx_t ctx)
{
    ca_csgn(res, x, GR_CA_CTX(ctx));

    if (ca_is_unknown(res, GR_CA_CTX(ctx)))
        return GR_UNABLE;

    return GR_SUCCESS;
}

int
_gr_ca_arg(ca_t res, const ca_t x, gr_ctx_t ctx)
{
    ca_arg(res, x, GR_CA_CTX(ctx));

    if (ca_is_unknown(res, GR_CA_CTX(ctx)))
        return GR_UNABLE;

    return GR_SUCCESS;
}


#define CMP_UNDEFINED -2
#define CMP_UNKNOWN -3

int
_ca_cmp(const ca_t x, const ca_t y, ca_ctx_t ctx);

int
_gr_ca_cmp(int * res, const ca_t x, const ca_t y, gr_ctx_t ctx)
{
    int cmp = _ca_cmp(x, y, GR_CA_CTX(ctx));

    if (cmp == CMP_UNKNOWN)
        return GR_UNABLE;

    if (cmp == CMP_UNDEFINED)
        return GR_DOMAIN;

    *res = cmp;
    return GR_SUCCESS;
}

int
_gr_ca_pi(ca_t res, gr_ctx_t ctx)
{
    if (ctx->which_ring == GR_CTX_REAL_ALGEBRAIC_CA ||
        ctx->which_ring == GR_CTX_COMPLEX_ALGEBRAIC_CA)
        return GR_DOMAIN;

    ca_pi(res, GR_CA_CTX(ctx));
    return handle_possible_special_value(res, ctx);
}

int
_gr_ca_exp(ca_t res, const ca_t x, gr_ctx_t ctx)
{
    if (ctx->which_ring == GR_CTX_REAL_ALGEBRAIC_CA ||
        ctx->which_ring == GR_CTX_COMPLEX_ALGEBRAIC_CA)
    {
        truth_t ok = ca_check_is_zero(x, GR_CA_CTX(ctx));

        if (ok == T_TRUE)
            return _gr_ca_one(res, ctx);

        return (ok == T_FALSE) ? GR_DOMAIN : GR_UNABLE;
    }
    else
    {
        ca_exp(res, x, GR_CA_CTX(ctx));
        return handle_possible_special_value(res, ctx);
    }
}

int
_gr_ca_sinh(ca_t res, const ca_t x, gr_ctx_t ctx)
{
    if (ctx->which_ring == GR_CTX_REAL_ALGEBRAIC_CA ||
        ctx->which_ring == GR_CTX_COMPLEX_ALGEBRAIC_CA)
    {
        truth_t ok = ca_check_is_zero(x, GR_CA_CTX(ctx));

        if (ok == T_TRUE)
            return _gr_ca_zero(res, ctx);

        return (ok == T_FALSE) ? GR_DOMAIN : GR_UNABLE;
    }
    else
    {
        int status = GR_SUCCESS;

        gr_ptr t, u;

        GR_TMP_INIT2(t, u, ctx);

        status |= gr_exp(t, x, ctx);
        status |= gr_inv(u, t, ctx);
        status |= gr_sub(res, t, u, ctx);
        status |= gr_mul_2exp_si(res, res, -1, ctx);

        GR_TMP_CLEAR2(t, u, ctx);

        return status | handle_possible_special_value(res, ctx);
    }
}

int
_gr_ca_cosh(ca_t res, const ca_t x, gr_ctx_t ctx)
{
    if (ctx->which_ring == GR_CTX_REAL_ALGEBRAIC_CA ||
        ctx->which_ring == GR_CTX_COMPLEX_ALGEBRAIC_CA)
    {
        truth_t ok = ca_check_is_zero(x, GR_CA_CTX(ctx));

        if (ok == T_TRUE)
            return _gr_ca_one(res, ctx);

        return (ok == T_FALSE) ? GR_DOMAIN : GR_UNABLE;
    }
    else
    {
        int status = GR_SUCCESS;

        gr_ptr t, u;

        GR_TMP_INIT2(t, u, ctx);

        status |= gr_exp(t, x, ctx);
        status |= gr_inv(u, t, ctx);
        status |= gr_add(res, t, u, ctx);
        status |= gr_mul_2exp_si(res, res, -1, ctx);

        GR_TMP_CLEAR2(t, u, ctx);

        return status | handle_possible_special_value(res, ctx);
    }
}

int
_gr_ca_tanh(ca_t res, const ca_t x, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;

    gr_ptr t, u;

    GR_TMP_INIT2(t, u, ctx);

    status |= gr_sinh(t, x, ctx);
    status |= gr_cosh(u, x, ctx);
    status |= gr_div(res, t, u, ctx);

    GR_TMP_CLEAR2(t, u, ctx);

    return status | handle_possible_special_value(res, ctx);
}

int
_gr_ca_coth(ca_t res, const ca_t x, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;

    gr_ptr t, u;

    GR_TMP_INIT2(t, u, ctx);

    status |= gr_sinh(t, x, ctx);
    status |= gr_cosh(u, x, ctx);
    status |= gr_div(res, u, t, ctx);

    GR_TMP_CLEAR2(t, u, ctx);

    return status | handle_possible_special_value(res, ctx);
}


int
_gr_ca_log(ca_t res, const ca_t x, gr_ctx_t ctx)
{
    if (ctx->which_ring == GR_CTX_REAL_ALGEBRAIC_CA ||
        ctx->which_ring == GR_CTX_COMPLEX_ALGEBRAIC_CA)
    {
        truth_t ok = ca_check_is_one(x, GR_CA_CTX(ctx));

        if (ok == T_TRUE)
            return _gr_ca_zero(res, ctx);

        return (ok == T_FALSE) ? GR_DOMAIN : GR_UNABLE;
    }
    else
    {
        ca_log(res, x, GR_CA_CTX(ctx));

        if (ctx->which_ring == GR_CTX_RR_CA)
        {
            truth_t ok = ca_check_is_real(res, GR_CA_CTX(ctx));

            if (ok == T_TRUE)
                return GR_SUCCESS;

            return (ok == T_FALSE) ? GR_DOMAIN : GR_UNABLE;
        }

        return handle_possible_special_value(res, ctx);
    }
}

int
_gr_ca_atan(ca_t res, const ca_t x, gr_ctx_t ctx)
{
    if (ctx->which_ring == GR_CTX_REAL_ALGEBRAIC_CA ||
        ctx->which_ring == GR_CTX_COMPLEX_ALGEBRAIC_CA)
    {
        truth_t ok = ca_check_is_zero(x, GR_CA_CTX(ctx));

        if (ok == T_TRUE)
            return _gr_ca_zero(res, ctx);

        return (ok == T_FALSE) ? GR_DOMAIN : GR_UNABLE;
    }
    else
    {
        ca_atan(res, x, GR_CA_CTX(ctx));
        return handle_possible_special_value(res, ctx);
    }
}

int
_gr_ca_sin(ca_t res, const ca_t x, gr_ctx_t ctx)
{
    if (ctx->which_ring == GR_CTX_REAL_ALGEBRAIC_CA ||
        ctx->which_ring == GR_CTX_COMPLEX_ALGEBRAIC_CA)
    {
        truth_t ok = ca_check_is_zero(x, GR_CA_CTX(ctx));

        if (ok == T_TRUE)
            return _gr_ca_zero(res, ctx);

        return (ok == T_FALSE) ? GR_DOMAIN : GR_UNABLE;
    }
    else
    {
        ca_sin(res, x, GR_CA_CTX(ctx));
        return handle_possible_special_value(res, ctx);
    }
}

int
_gr_ca_cos(ca_t res, const ca_t x, gr_ctx_t ctx)
{
    if (ctx->which_ring == GR_CTX_REAL_ALGEBRAIC_CA ||
        ctx->which_ring == GR_CTX_COMPLEX_ALGEBRAIC_CA)
    {
        truth_t ok = ca_check_is_zero(x, GR_CA_CTX(ctx));

        if (ok == T_TRUE)
            return _gr_ca_one(res, ctx);

        return (ok == T_FALSE) ? GR_DOMAIN : GR_UNABLE;
    }
    else
    {
        ca_cos(res, x, GR_CA_CTX(ctx));
        return handle_possible_special_value(res, ctx);
    }
}

int
_gr_ca_tan(ca_t res, const ca_t x, gr_ctx_t ctx)
{
    if (ctx->which_ring == GR_CTX_REAL_ALGEBRAIC_CA ||
        ctx->which_ring == GR_CTX_COMPLEX_ALGEBRAIC_CA)
    {
        truth_t ok = ca_check_is_zero(x, GR_CA_CTX(ctx));

        if (ok == T_TRUE)
            return _gr_ca_zero(res, ctx);

        return (ok == T_FALSE) ? GR_DOMAIN : GR_UNABLE;
    }
    else
    {
        ca_tan(res, x, GR_CA_CTX(ctx));
        return handle_possible_special_value(res, ctx);
    }
}

int
_gr_ca_cot(ca_t res, const ca_t x, gr_ctx_t ctx)
{
    if (ctx->which_ring == GR_CTX_REAL_ALGEBRAIC_CA ||
        ctx->which_ring == GR_CTX_COMPLEX_ALGEBRAIC_CA)
    {
        return GR_DOMAIN;
    }
    else
    {
        ca_cot(res, x, GR_CA_CTX(ctx));
        return handle_possible_special_value(res, ctx);
    }
}

int
_gr_ca_asin(ca_t res, const ca_t x, gr_ctx_t ctx)
{
    if (ctx->which_ring == GR_CTX_REAL_ALGEBRAIC_CA ||
        ctx->which_ring == GR_CTX_COMPLEX_ALGEBRAIC_CA)
    {
        truth_t ok = ca_check_is_zero(x, GR_CA_CTX(ctx));

        if (ok == T_TRUE)
            return _gr_ca_zero(res, ctx);

        return (ok == T_FALSE) ? GR_DOMAIN : GR_UNABLE;
    }
    else
    {
        ca_asin(res, x, GR_CA_CTX(ctx));
        return handle_possible_special_value(res, ctx);
    }
}

int
_gr_ca_acos(ca_t res, const ca_t x, gr_ctx_t ctx)
{
    if (ctx->which_ring == GR_CTX_REAL_ALGEBRAIC_CA ||
        ctx->which_ring == GR_CTX_COMPLEX_ALGEBRAIC_CA)
    {
        truth_t ok = ca_check_is_one(x, GR_CA_CTX(ctx));

        if (ok == T_TRUE)
            return _gr_ca_zero(res, ctx);

        return (ok == T_FALSE) ? GR_DOMAIN : GR_UNABLE;
    }
    else
    {
        ca_acos(res, x, GR_CA_CTX(ctx));
        return handle_possible_special_value(res, ctx);
    }
}

int
_gr_ca_erf(ca_t res, const ca_t x, gr_ctx_t ctx)
{
    if (ctx->which_ring == GR_CTX_REAL_ALGEBRAIC_CA ||
        ctx->which_ring == GR_CTX_COMPLEX_ALGEBRAIC_CA)
    {
        truth_t ok = ca_check_is_zero(x, GR_CA_CTX(ctx));

        if (ok == T_TRUE)
            return _gr_ca_zero(res, ctx);

        return GR_UNABLE;
    }
    else
    {
        ca_erf(res, x, GR_CA_CTX(ctx));
        return handle_possible_special_value(res, ctx);
    }
}

int
_gr_ca_erfi(ca_t res, const ca_t x, gr_ctx_t ctx)
{
    if (ctx->which_ring == GR_CTX_REAL_ALGEBRAIC_CA ||
        ctx->which_ring == GR_CTX_COMPLEX_ALGEBRAIC_CA)
    {
        truth_t ok = ca_check_is_zero(x, GR_CA_CTX(ctx));

        if (ok == T_TRUE)
            return _gr_ca_zero(res, ctx);

        return GR_UNABLE;
    }
    else
    {
        ca_erfi(res, x, GR_CA_CTX(ctx));
        return handle_possible_special_value(res, ctx);
    }
}

int
_gr_ca_erfc(ca_t res, const ca_t x, gr_ctx_t ctx)
{
    if (ctx->which_ring == GR_CTX_REAL_ALGEBRAIC_CA ||
        ctx->which_ring == GR_CTX_COMPLEX_ALGEBRAIC_CA)
    {
        truth_t ok = ca_check_is_zero(x, GR_CA_CTX(ctx));

        if (ok == T_TRUE)
            return _gr_ca_one(res, ctx);

        return GR_UNABLE;
    }
    else
    {
        ca_erfc(res, x, GR_CA_CTX(ctx));
        return handle_possible_special_value(res, ctx);
    }
}

int
_gr_ca_gamma(ca_t res, const ca_t x, gr_ctx_t ctx)
{
    if (ctx->which_ring == GR_CTX_REAL_ALGEBRAIC_CA ||
        ctx->which_ring == GR_CTX_COMPLEX_ALGEBRAIC_CA)
    {
        truth_t ok = ca_check_is_integer(x, GR_CA_CTX(ctx));

        if (ok != T_TRUE)
            return GR_UNABLE;
    }

    {
        ca_gamma(res, x, GR_CA_CTX(ctx));
        return handle_possible_special_value(res, ctx);
    }
}

int
_gr_ca_poly_mullow(ca_ptr res,
    ca_srcptr poly1, slong len1,
    ca_srcptr poly2, slong len2, slong n, gr_ctx_t ctx)
{
    _ca_poly_mullow(res, poly1, len1, poly2, len2, n, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_poly_roots(gr_vec_t roots, gr_vec_t mult, const gr_poly_t poly, int flags, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    ca_vec_t ca_roots;
    ulong * exp;
    slong deg;
    gr_ctx_t ZZ;
    slong i;
    fmpz_t exp_z;

    if (poly->length == 0)
        return GR_DOMAIN;

    deg = poly->length - 1;

    gr_ctx_init_fmpz(ZZ);
    fmpz_init(exp_z);

    ca_vec_init(ca_roots, 0, GR_CA_CTX(ctx));
    exp = flint_malloc(sizeof(ulong) * deg);

    if (ca_poly_roots(ca_roots, exp, (const ca_poly_struct *) poly, GR_CA_CTX(ctx)))
    {
        gr_vec_set_length(roots, 0, ctx);
        gr_vec_set_length(mult, 0, ZZ);

        for (i = 0; i < ca_roots->length; i++)
        {
            if (ctx->which_ring == GR_CTX_REAL_ALGEBRAIC_CA || ctx->which_ring == GR_CTX_RR_CA)
            {
                truth_t is_real = ca_check_is_real(ca_roots->entries + i, GR_CA_CTX(ctx));

                if (is_real == T_FALSE)
                    continue;

                if (is_real == T_UNKNOWN)
                {
                    status = GR_UNABLE;
                    break;
                }
            }

            fmpz_set_ui(exp_z, exp[i]);
            status |= gr_vec_append(roots, ca_roots->entries + i, ctx);
            status |= gr_vec_append(mult, exp_z, ZZ);
        }
    }
    else
    {
        status = GR_UNABLE;
        gr_vec_set_length(roots, 0, ctx);
        gr_vec_set_length(mult, 0, ZZ);
    }

    ca_vec_clear(ca_roots, GR_CA_CTX(ctx));
    flint_free(exp);

    gr_ctx_clear(ZZ);
    fmpz_clear(exp_z);

/*
    if (status == GR_SUCCESS)
    {
        _acb_vec_sort_pretty(croots, deg);

        if (arb_roots)
        {

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
*/

    return status;
}

int
_gr_ca_mat_mul(ca_mat_t res, const ca_mat_t x, const ca_mat_t y, gr_ctx_t ctx)
{
    ca_mat_mul(res, x, y, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_mat_det(ca_t res, const ca_mat_t x, gr_ctx_t ctx)
{
    ca_mat_det(res, x, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_ctx_clear(gr_ctx_t ctx)
{
    ca_ctx_clear(GR_CA_CTX(ctx));
    flint_free(GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

truth_t
_gr_ca_ctx_is_field(gr_ctx_t ctx)
{
    return (ctx->which_ring == GR_CTX_COMPLEX_EXTENDED_CA) ? T_FALSE : T_TRUE;
}


truth_t
_gr_ca_ctx_is_algebraically_closed(gr_ctx_t ctx)
{
    return (ctx->which_ring == GR_CTX_CC_CA || ctx->which_ring == GR_CTX_COMPLEX_ALGEBRAIC_CA) ? T_TRUE : T_FALSE;
}

truth_t
_gr_ca_ctx_is_ordered_ring(gr_ctx_t ctx)
{
    return (ctx->which_ring == GR_CTX_RR_CA || ctx->which_ring == GR_CTX_REAL_ALGEBRAIC_CA) ? T_TRUE : T_FALSE;
}

int
_gr_ca_mat_find_nonzero_pivot(slong * pivot_row, ca_mat_t mat, slong start_row, slong end_row, slong column, gr_ctx_t ctx)
{
    truth_t ok;

    ok = ca_mat_find_pivot(pivot_row, mat, start_row, end_row, column, GR_CA_CTX(ctx));

    if (ok == T_TRUE)
        return GR_SUCCESS;
    else if (ok == T_FALSE)
        return GR_DOMAIN;
    else
        return GR_UNABLE;
}

int _ca_methods_initialized = 0;

FLINT_DLL gr_static_method_table _ca_methods;

gr_method_tab_input _ca_methods_input[] =
{
    {GR_METHOD_CTX_CLEAR,       (gr_funcptr) _gr_ca_ctx_clear},
    {GR_METHOD_CTX_WRITE,       (gr_funcptr) _gr_ca_ctx_write},

    {GR_METHOD_CTX_IS_RING,     (gr_funcptr) _gr_ca_ctx_is_field},
    {GR_METHOD_CTX_IS_COMMUTATIVE_RING, (gr_funcptr) _gr_ca_ctx_is_field},
    {GR_METHOD_CTX_IS_INTEGRAL_DOMAIN,  (gr_funcptr) _gr_ca_ctx_is_field},
    {GR_METHOD_CTX_IS_FIELD,            (gr_funcptr) _gr_ca_ctx_is_field},

    {GR_METHOD_CTX_IS_UNIQUE_FACTORIZATION_DOMAIN, (gr_funcptr) _gr_ca_ctx_is_field},
    {GR_METHOD_CTX_IS_FINITE, (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_FINITE_CHARACTERISTIC, (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_ALGEBRAICALLY_CLOSED, (gr_funcptr) _gr_ca_ctx_is_algebraically_closed},
    {GR_METHOD_CTX_IS_ORDERED_RING, (gr_funcptr) _gr_ca_ctx_is_ordered_ring},

    /* important */
    {GR_METHOD_CTX_IS_THREADSAFE,       (gr_funcptr) gr_generic_ctx_predicate_false},

    {GR_METHOD_CTX_IS_EXACT,     (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_CANONICAL, (gr_funcptr) gr_generic_ctx_predicate_false},

    {GR_METHOD_INIT,            (gr_funcptr) _gr_ca_init},

    {GR_METHOD_CLEAR,           (gr_funcptr) _gr_ca_clear},
    {GR_METHOD_SWAP,            (gr_funcptr) _gr_ca_swap},
    {GR_METHOD_SET_SHALLOW,     (gr_funcptr) _gr_ca_set_shallow},
    {GR_METHOD_RANDTEST,        (gr_funcptr) _gr_ca_randtest},
    {GR_METHOD_WRITE,           (gr_funcptr) _gr_ca_write},

    {GR_METHOD_ZERO,            (gr_funcptr) _gr_ca_zero},
    {GR_METHOD_ONE,             (gr_funcptr) _gr_ca_one},
    {GR_METHOD_IS_ZERO,         (gr_funcptr) _gr_ca_is_zero},
    {GR_METHOD_IS_ONE,          (gr_funcptr) _gr_ca_is_one},
    {GR_METHOD_IS_NEG_ONE,      (gr_funcptr) _gr_ca_is_neg_one},
    {GR_METHOD_EQUAL,           (gr_funcptr) _gr_ca_equal},

    {GR_METHOD_SET,             (gr_funcptr) _gr_ca_set},
    {GR_METHOD_SET_SI,          (gr_funcptr) _gr_ca_set_si},
    {GR_METHOD_SET_UI,          (gr_funcptr) _gr_ca_set_ui},
    {GR_METHOD_SET_FMPZ,        (gr_funcptr) _gr_ca_set_fmpz},
    {GR_METHOD_SET_FMPQ,        (gr_funcptr) _gr_ca_set_fmpq},
    {GR_METHOD_SET_OTHER,       (gr_funcptr) _gr_ca_set_other},

    {GR_METHOD_GET_SI,          (gr_funcptr) _gr_ca_get_si},
    {GR_METHOD_GET_UI,          (gr_funcptr) _gr_ca_get_ui},
    {GR_METHOD_GET_FMPZ,        (gr_funcptr) _gr_ca_get_fmpz},
    {GR_METHOD_GET_FMPQ,        (gr_funcptr) _gr_ca_get_fmpq},
    {GR_METHOD_GET_D,           (gr_funcptr) _gr_ca_get_d},

    {GR_METHOD_GET_FEXPR,                 (gr_funcptr) _gr_ca_get_fexpr},
    {GR_METHOD_GET_FEXPR_SERIALIZE,       (gr_funcptr) _gr_ca_get_fexpr_serialize},

    {GR_METHOD_NEG,             (gr_funcptr) _gr_ca_neg},

    {GR_METHOD_ADD,             (gr_funcptr) _gr_ca_add},
    {GR_METHOD_ADD_UI,          (gr_funcptr) _gr_ca_add_ui},
    {GR_METHOD_ADD_SI,          (gr_funcptr) _gr_ca_add_si},
    {GR_METHOD_ADD_FMPZ,        (gr_funcptr) _gr_ca_add_fmpz},
    {GR_METHOD_ADD_FMPQ,        (gr_funcptr) _gr_ca_add_fmpq},

    {GR_METHOD_SUB,             (gr_funcptr) _gr_ca_sub},
    {GR_METHOD_SUB_UI,          (gr_funcptr) _gr_ca_sub_ui},
    {GR_METHOD_SUB_SI,          (gr_funcptr) _gr_ca_sub_si},
    {GR_METHOD_SUB_FMPZ,        (gr_funcptr) _gr_ca_sub_fmpz},
    {GR_METHOD_SUB_FMPQ,        (gr_funcptr) _gr_ca_sub_fmpq},

    {GR_METHOD_MUL,             (gr_funcptr) _gr_ca_mul},
    {GR_METHOD_MUL_UI,          (gr_funcptr) _gr_ca_mul_ui},
    {GR_METHOD_MUL_SI,          (gr_funcptr) _gr_ca_mul_si},
    {GR_METHOD_MUL_FMPZ,        (gr_funcptr) _gr_ca_mul_fmpz},
    {GR_METHOD_MUL_FMPQ,        (gr_funcptr) _gr_ca_mul_fmpq},

    {GR_METHOD_DIV,             (gr_funcptr) _gr_ca_div},
    {GR_METHOD_DIV_UI,          (gr_funcptr) _gr_ca_div_ui},
    {GR_METHOD_DIV_SI,          (gr_funcptr) _gr_ca_div_si},
    {GR_METHOD_DIV_FMPZ,        (gr_funcptr) _gr_ca_div_fmpz},
    {GR_METHOD_DIV_FMPQ,        (gr_funcptr) _gr_ca_div_fmpq},

    {GR_METHOD_IS_INVERTIBLE,   (gr_funcptr) _gr_ca_is_invertible},
    {GR_METHOD_INV,             (gr_funcptr) _gr_ca_inv},

    {GR_METHOD_POW,             (gr_funcptr) _gr_ca_pow},
    {GR_METHOD_POW_UI,          (gr_funcptr) _gr_ca_pow_ui},
    {GR_METHOD_POW_SI,          (gr_funcptr) _gr_ca_pow_si},
    {GR_METHOD_POW_FMPZ,        (gr_funcptr) _gr_ca_pow_fmpz},
    {GR_METHOD_POW_FMPQ,        (gr_funcptr) _gr_ca_pow_fmpq},

    {GR_METHOD_IS_SQUARE,       (gr_funcptr) _gr_ca_is_square},
    {GR_METHOD_SQRT,            (gr_funcptr) _gr_ca_sqrt},
    {GR_METHOD_RSQRT,           (gr_funcptr) _gr_ca_rsqrt},

    {GR_METHOD_POS_INF,         (gr_funcptr) _gr_ca_pos_inf},
    {GR_METHOD_NEG_INF,         (gr_funcptr) _gr_ca_neg_inf},
    {GR_METHOD_UINF,            (gr_funcptr) _gr_ca_uinf},
    {GR_METHOD_UNDEFINED,       (gr_funcptr) _gr_ca_undefined},
    {GR_METHOD_UNKNOWN,         (gr_funcptr) _gr_ca_unknown},

    {GR_METHOD_FLOOR,           (gr_funcptr) _gr_ca_floor},
    {GR_METHOD_CEIL,            (gr_funcptr) _gr_ca_ceil},
    {GR_METHOD_TRUNC,           (gr_funcptr) _gr_ca_trunc},
    {GR_METHOD_NINT,            (gr_funcptr) _gr_ca_nint},

    {GR_METHOD_I,               (gr_funcptr) _gr_ca_i},
    {GR_METHOD_ABS,             (gr_funcptr) _gr_ca_abs},
    {GR_METHOD_CONJ,            (gr_funcptr) _gr_ca_conj},
    {GR_METHOD_RE,              (gr_funcptr) _gr_ca_re},
    {GR_METHOD_IM,              (gr_funcptr) _gr_ca_im},
    {GR_METHOD_SGN,             (gr_funcptr) _gr_ca_sgn},
    {GR_METHOD_CSGN,            (gr_funcptr) _gr_ca_csgn},
    {GR_METHOD_ARG,             (gr_funcptr) _gr_ca_arg},

    {GR_METHOD_CMP,              (gr_funcptr) _gr_ca_cmp},

    {GR_METHOD_PI,              (gr_funcptr) _gr_ca_pi},
    {GR_METHOD_EXP,             (gr_funcptr) _gr_ca_exp},
    {GR_METHOD_LOG,             (gr_funcptr) _gr_ca_log},

    {GR_METHOD_SIN,             (gr_funcptr) _gr_ca_sin},
    {GR_METHOD_COS,             (gr_funcptr) _gr_ca_cos},
    {GR_METHOD_TAN,             (gr_funcptr) _gr_ca_tan},
    {GR_METHOD_COT,             (gr_funcptr) _gr_ca_cot},

    {GR_METHOD_SINH,             (gr_funcptr) _gr_ca_sinh},
    {GR_METHOD_COSH,             (gr_funcptr) _gr_ca_cosh},
    {GR_METHOD_TANH,             (gr_funcptr) _gr_ca_tanh},
    {GR_METHOD_COTH,             (gr_funcptr) _gr_ca_coth},

    {GR_METHOD_ATAN,            (gr_funcptr) _gr_ca_atan},
    {GR_METHOD_ASIN,            (gr_funcptr) _gr_ca_asin},
    {GR_METHOD_ACOS,            (gr_funcptr) _gr_ca_acos},

    {GR_METHOD_ERF,             (gr_funcptr) _gr_ca_erf},
    {GR_METHOD_ERFC,            (gr_funcptr) _gr_ca_erfc},
    {GR_METHOD_ERFI,            (gr_funcptr) _gr_ca_erfi},

    {GR_METHOD_GAMMA,            (gr_funcptr) _gr_ca_gamma},

    {GR_METHOD_POLY_MULLOW,     (gr_funcptr) _gr_ca_poly_mullow},
    {GR_METHOD_POLY_ROOTS,      (gr_funcptr) _gr_ca_poly_roots},
    /* {GR_METHOD_POLY_ROOTS_OTHER,(gr_funcptr) _gr_ca_poly_roots_other}, */
    {GR_METHOD_MAT_MUL,         (gr_funcptr) _gr_ca_mat_mul},
    {GR_METHOD_MAT_DET,         (gr_funcptr) _gr_ca_mat_det},
    {GR_METHOD_MAT_FIND_NONZERO_PIVOT,      (gr_funcptr) _gr_ca_mat_find_nonzero_pivot},

    {0,                         (gr_funcptr) NULL},
};

void
_gr_ctx_init_ca(gr_ctx_t ctx, int which_ring)
{
    ctx->which_ring = which_ring;
    ctx->sizeof_elem = sizeof(ca_struct);
    ctx->size_limit = WORD_MAX;

    GR_CTX_DATA_AS_PTR(ctx) = flint_malloc(sizeof(ca_ctx_struct));
    ca_ctx_init(GR_CA_CTX(ctx));

    ctx->methods = _ca_methods;

    if (!_ca_methods_initialized)
    {
        gr_method_tab_init(_ca_methods, _ca_methods_input);
        _ca_methods_initialized = 1;
    }
}

/* ca_ctx is a (ca_ctx_struct *) */
void
_gr_ctx_init_ca_from_ref(gr_ctx_t ctx, int which_ring, void * ca_ctx)
{
    ctx->which_ring = which_ring;
    ctx->sizeof_elem = sizeof(ca_struct);
    ctx->size_limit = WORD_MAX;

    GR_CTX_DATA_AS_PTR(ctx) = ca_ctx;

    ctx->methods = _ca_methods;

    if (!_ca_methods_initialized)
    {
        gr_method_tab_init(_ca_methods, _ca_methods_input);
        _ca_methods_initialized = 1;
    }
}

void
gr_ctx_init_real_ca(gr_ctx_t ctx)
{
    _gr_ctx_init_ca(ctx, GR_CTX_RR_CA);
}

void
gr_ctx_init_complex_ca(gr_ctx_t ctx)
{
    _gr_ctx_init_ca(ctx, GR_CTX_CC_CA);
}

void
gr_ctx_init_real_algebraic_ca(gr_ctx_t ctx)
{
    _gr_ctx_init_ca(ctx, GR_CTX_REAL_ALGEBRAIC_CA);
}

void
gr_ctx_init_complex_algebraic_ca(gr_ctx_t ctx)
{
    _gr_ctx_init_ca(ctx, GR_CTX_COMPLEX_ALGEBRAIC_CA);
}

void
gr_ctx_init_complex_extended_ca(gr_ctx_t ctx)
{
    _gr_ctx_init_ca(ctx, GR_CTX_COMPLEX_EXTENDED_CA);
}
