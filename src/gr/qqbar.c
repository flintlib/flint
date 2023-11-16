/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"
#include "fmpz_poly_factor.h"
#include "fmpq.h"
#include "fmpzi.h"
#include "fexpr.h"
#include "qqbar.h"
#include "gr.h"
#include "gr_generic.h"
#include "gr_vec.h"
#include "gr_poly.h"
#include "ca.h"

typedef struct
{
    /* todo: use which_ring */
    int real_only;      /* field restricted to real algebraics instead of complex? */
    slong deg_limit;    /* todo */
    slong bits_limit;   /* todo */
}
gr_qqbar_ctx;

#define QQBAR_CTX(ring_ctx) ((gr_qqbar_ctx *)(ring_ctx))

int
_gr_qqbar_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    if (QQBAR_CTX(ctx)->real_only)
        gr_stream_write(out, "Real algebraic numbers (qqbar)");
    else
        gr_stream_write(out, "Complex algebraic numbers (qqbar)");
    return GR_SUCCESS;
}

void
_gr_qqbar_init(qqbar_t x, const gr_ctx_t ctx)
{
    qqbar_init(x);
}

void
_gr_qqbar_clear(qqbar_t x, const gr_ctx_t ctx)
{
    qqbar_clear(x);
}

void
_gr_qqbar_swap(qqbar_t x, qqbar_t y, const gr_ctx_t ctx)
{
    qqbar_t t;
    *t = *x;
    *x = *y;
    *y = *t;
}

void
_gr_qqbar_set_shallow(qqbar_t res, const qqbar_t x, const gr_ctx_t ctx)
{
    *res = *x;
}

/* todo: limits */
int
_gr_qqbar_randtest(qqbar_t res, flint_rand_t state, const gr_ctx_t ctx)
{
    slong deg_limit, bits_limit;
    int rcase;

    rcase = n_randint(state, 10);

    if (rcase == 0)
    {
        deg_limit = 4;
        bits_limit = 10;
    }
    else if (rcase <= 3)
    {
        deg_limit = 2;
        bits_limit = 10;
    }
    else
    {
        deg_limit = 1;
        bits_limit = 10;
    }

    if (QQBAR_CTX(ctx)->real_only)
        qqbar_randtest_real(res, state, deg_limit, bits_limit);
    else
        qqbar_randtest(res, state, deg_limit, bits_limit);

    return GR_SUCCESS;
}

/* todo: different styles */


void
qqbar_get_decimal_root_nearest(char ** re_s, char ** im_s, const qqbar_t x, slong default_digits);

int
_gr_qqbar_write(gr_stream_t out, const qqbar_t x, const gr_ctx_t ctx)
{
    char *re_s, *im_s;

    if (qqbar_is_rational(x))
    {
        fmpq_t t;
        fmpq_init(t);
        qqbar_get_fmpq(t, x);
        gr_stream_write_fmpz(out, fmpq_numref(t));
        if (!fmpz_is_one(fmpq_denref(t)))
        {
            gr_stream_write(out, "/");
            gr_stream_write_fmpz(out, fmpq_denref(t));
        }
        fmpq_clear(t);
    }
    else
    {
        qqbar_get_decimal_root_nearest(&re_s, &im_s, x, 6);

        gr_stream_write(out, "Root a = ");

        if (re_s != NULL)
        {
            gr_stream_write_free(out, re_s);
        }

        if (im_s != NULL)
        {
            if (re_s != NULL)
            {
                if (im_s[0] == '-')
                {
                    gr_stream_write(out, " - ");
                    gr_stream_write(out, im_s + 1);
                    flint_free(im_s);
                }
                else
                {
                    gr_stream_write(out, " + ");
                    gr_stream_write_free(out, im_s);
                }
            }
            else
            {
                gr_stream_write_free(out, im_s);
            }

            gr_stream_write(out, "*I");
        }
        gr_stream_write(out, " of ");
        gr_stream_write_free(out, fmpz_poly_get_str_pretty(QQBAR_POLY(x), "a"));
    }

    return GR_SUCCESS;
}

int
_gr_qqbar_zero(qqbar_t x, const gr_ctx_t ctx)
{
    qqbar_zero(x);
    return GR_SUCCESS;
}

int
_gr_qqbar_one(qqbar_t x, const gr_ctx_t ctx)
{
    qqbar_one(x);
    return GR_SUCCESS;
}

int
_gr_qqbar_set_si(qqbar_t res, slong v, const gr_ctx_t ctx)
{
    qqbar_set_si(res, v);
    return GR_SUCCESS;
}

int
_gr_qqbar_set_ui(qqbar_t res, ulong v, const gr_ctx_t ctx)
{
    qqbar_set_ui(res, v);
    return GR_SUCCESS;
}

int
_gr_qqbar_set_fmpz(qqbar_t res, const fmpz_t v, const gr_ctx_t ctx)
{
    qqbar_set_fmpz(res, v);
    return GR_SUCCESS;
}

int
_gr_qqbar_set_fmpq(qqbar_t res, const fmpq_t v, const gr_ctx_t ctx)
{
    qqbar_set_fmpq(res, v);
    return GR_SUCCESS;
}

int
_gr_qqbar_set_d(qqbar_t res, double v, const gr_ctx_t ctx)
{
    if (qqbar_set_d(res, v))
        return GR_SUCCESS;
    else
        return GR_DOMAIN;
}

truth_t
_gr_qqbar_is_zero(const qqbar_t x, const gr_ctx_t ctx)
{
    return qqbar_is_zero(x) ? T_TRUE : T_FALSE;
}

truth_t
_gr_qqbar_is_one(const qqbar_t x, const gr_ctx_t ctx)
{
    return qqbar_is_one(x) ? T_TRUE : T_FALSE;
}

truth_t
_gr_qqbar_is_neg_one(const qqbar_t x, const gr_ctx_t ctx)
{
    return qqbar_is_neg_one(x) ? T_TRUE : T_FALSE;
}

truth_t
_gr_qqbar_equal(const qqbar_t x, const qqbar_t y, const gr_ctx_t ctx)
{
    return qqbar_equal(x, y) ? T_TRUE : T_FALSE;
}

int
_gr_qqbar_set(qqbar_t res, const qqbar_t x, const gr_ctx_t ctx)
{
    qqbar_set(res, x);
    return GR_SUCCESS;
}

/* todo: move */
void
qqbar_set_fmpzi(qqbar_t res, const fmpzi_t x)
{
    if (fmpz_is_zero(fmpzi_imagref(x)))
    {
        qqbar_set_fmpz(res, fmpzi_realref(x));
    }
    else
    {
        fmpz_poly_fit_length(QQBAR_POLY(res), 3);
        _fmpz_poly_set_length(QQBAR_POLY(res), 3);
        fmpzi_norm(QQBAR_COEFFS(res), x);
        fmpz_mul_si(QQBAR_COEFFS(res) + 1, fmpzi_realref(x), -2);
        fmpz_one(QQBAR_COEFFS(res) + 2);
        arb_set_round_fmpz(acb_realref(QQBAR_ENCLOSURE(res)), fmpzi_realref(x), QQBAR_DEFAULT_PREC);
        arb_set_round_fmpz(acb_imagref(QQBAR_ENCLOSURE(res)), fmpzi_imagref(x), QQBAR_DEFAULT_PREC);
    }
}

int
_gr_qqbar_set_other(qqbar_t res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
{
    switch (x_ctx->which_ring)
    {
        case GR_CTX_FMPZ:
            qqbar_set_fmpz(res, x);
            return GR_SUCCESS;

        case GR_CTX_FMPQ:
            qqbar_set_fmpq(res, x);
            return GR_SUCCESS;

        case GR_CTX_FMPZI:
            if (ctx->which_ring == GR_CTX_REAL_ALGEBRAIC_QQBAR && !fmpz_is_zero(fmpzi_imagref((const fmpzi_struct *) x)))
                return GR_DOMAIN;
            qqbar_set_fmpzi(res, x);
            return GR_SUCCESS;

        case GR_CTX_REAL_ALGEBRAIC_QQBAR:
        case GR_CTX_COMPLEX_ALGEBRAIC_QQBAR:
            if (ctx->which_ring == GR_CTX_REAL_ALGEBRAIC_QQBAR && !qqbar_is_real(x))
                return GR_DOMAIN;
            qqbar_set(res, x);
            return GR_SUCCESS;

        /* todo: all cases */
        case GR_CTX_RR_CA:
        case GR_CTX_CC_CA:
        case GR_CTX_REAL_ALGEBRAIC_CA:
        case GR_CTX_COMPLEX_ALGEBRAIC_CA:
        case GR_CTX_COMPLEX_EXTENDED_CA:

            if (!ca_get_qqbar(res, x, gr_ctx_data_as_ptr(x_ctx)))
                return GR_UNABLE;

            if (ctx->which_ring == GR_CTX_COMPLEX_ALGEBRAIC_QQBAR ||
                x_ctx->which_ring == GR_CTX_REAL_ALGEBRAIC_CA ||
                x_ctx->which_ring == GR_CTX_RR_CA)
                return GR_SUCCESS;

            if (!qqbar_is_real(res))
            {
                qqbar_zero(res);
                return GR_DOMAIN;
            }

            return GR_SUCCESS;

    }

    return gr_generic_set_other(res, x, x_ctx, ctx);
}

int
_gr_qqbar_get_fmpz(fmpz_t res, const qqbar_t x, const gr_ctx_t ctx)
{
    if (!qqbar_is_integer(x))
        return GR_DOMAIN;

    qqbar_get_fmpz(res, x);
    return GR_SUCCESS;
}

int
_gr_qqbar_get_ui(ulong * res, const qqbar_t x, const gr_ctx_t ctx)
{
    fmpz_t t;
    int status;

    if (!qqbar_is_integer(x))
        return GR_DOMAIN;

    fmpz_init(t);
    qqbar_get_fmpz(t, x);

    if (fmpz_sgn(t) < 0 || fmpz_cmp_ui(t, UWORD_MAX) > 0)
    {
        status = GR_DOMAIN;
    }
    else
    {
        *res = fmpz_get_ui(t);
        status = GR_SUCCESS;
    }

    fmpz_clear(t);
    return status;
}

int
_gr_qqbar_get_si(slong * res, const qqbar_t x, const gr_ctx_t ctx)
{
    fmpz_t t;
    int status;

    if (!qqbar_is_integer(x))
        return GR_DOMAIN;

    fmpz_init(t);
    qqbar_get_fmpz(t, x);

    if (!fmpz_fits_si(t))
    {
        status = GR_DOMAIN;
    }
    else
    {
        *res = fmpz_get_si(t);
        status = GR_SUCCESS;
    }

    fmpz_clear(t);
    return status;
}

int
_gr_qqbar_get_fmpq(fmpq_t res, const qqbar_t x, const gr_ctx_t ctx)
{
    if (!qqbar_is_rational(x))
        return GR_DOMAIN;

    qqbar_get_fmpq(res, x);
    return GR_SUCCESS;
}

int
_gr_qqbar_get_d(double * res, const qqbar_t x, const gr_ctx_t ctx)
{
    arb_t t;

    if (!qqbar_is_real(x))
        return GR_DOMAIN;

    arb_init(t);
    qqbar_get_arb(t, x, 64);
    *res = arf_get_d(arb_midref(t), ARF_RND_NEAR);
    arb_clear(t);
    return GR_SUCCESS;
}

int
_gr_qqbar_get_fexpr(fexpr_t res, const qqbar_t x, const gr_ctx_t ctx)
{
    if (!qqbar_get_fexpr_formula(res, x, QQBAR_FORMULA_GAUSSIANS | QQBAR_FORMULA_QUADRATICS))
        qqbar_get_fexpr_root_nearest(res, x);
    return GR_SUCCESS;
}

int
_gr_qqbar_get_fexpr_serialize(fexpr_t res, const qqbar_t x, const gr_ctx_t ctx)
{
    qqbar_get_fexpr_repr(res, x);
    return GR_SUCCESS;
}


/* todo */
int _gr_qqbar_set_fexpr(gr_ptr res, fexpr_vec_t inputs, gr_vec_t outputs, const fexpr_t x, gr_ctx_t ctx)
{
    if (qqbar_set_fexpr(res, x))
    {
        if (QQBAR_CTX(ctx)->real_only && !qqbar_is_real(res))
            return GR_DOMAIN;
        else
            return GR_SUCCESS;
    }
    else
    {
        return gr_generic_set_fexpr(res, inputs, outputs, x, ctx);
    }
}

int
_gr_qqbar_neg(qqbar_t res, const qqbar_t x, const gr_ctx_t ctx)
{
    qqbar_neg(res, x);
    return GR_SUCCESS;
}

int
_gr_qqbar_add(qqbar_t res, const qqbar_t x, const qqbar_t y, const gr_ctx_t ctx)
{
    qqbar_add(res, x, y);
    return GR_SUCCESS;
}

int
_gr_qqbar_add_si(qqbar_t res, const qqbar_t x, slong y, const gr_ctx_t ctx)
{
    qqbar_add_si(res, x, y);
    return GR_SUCCESS;
}

int
_gr_qqbar_add_ui(qqbar_t res, const qqbar_t x, ulong y, const gr_ctx_t ctx)
{
    qqbar_add_ui(res, x, y);
    return GR_SUCCESS;
}

int
_gr_qqbar_add_fmpz(qqbar_t res, const qqbar_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    qqbar_add_fmpz(res, x, y);
    return GR_SUCCESS;
}

int
_gr_qqbar_add_fmpq(qqbar_t res, const qqbar_t x, const fmpq_t y, const gr_ctx_t ctx)
{
    qqbar_add_fmpq(res, x, y);
    return GR_SUCCESS;
}

int
_gr_qqbar_sub(qqbar_t res, const qqbar_t x, const qqbar_t y, const gr_ctx_t ctx)
{
    qqbar_sub(res, x, y);
    return GR_SUCCESS;
}

int
_gr_qqbar_sub_si(qqbar_t res, const qqbar_t x, slong y, const gr_ctx_t ctx)
{
    qqbar_sub_si(res, x, y);
    return GR_SUCCESS;
}

int
_gr_qqbar_sub_ui(qqbar_t res, const qqbar_t x, ulong y, const gr_ctx_t ctx)
{
    qqbar_sub_ui(res, x, y);
    return GR_SUCCESS;
}

int
_gr_qqbar_sub_fmpz(qqbar_t res, const qqbar_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    qqbar_sub_fmpz(res, x, y);
    return GR_SUCCESS;
}

int
_gr_qqbar_sub_fmpq(qqbar_t res, const qqbar_t x, const fmpq_t y, const gr_ctx_t ctx)
{
    qqbar_sub_fmpq(res, x, y);
    return GR_SUCCESS;
}

int
_gr_qqbar_mul(qqbar_t res, const qqbar_t x, const qqbar_t y, const gr_ctx_t ctx)
{
    qqbar_mul(res, x, y);
    return GR_SUCCESS;
}

int
_gr_qqbar_mul_si(qqbar_t res, const qqbar_t x, slong y, const gr_ctx_t ctx)
{
    qqbar_mul_si(res, x, y);
    return GR_SUCCESS;
}

int
_gr_qqbar_mul_ui(qqbar_t res, const qqbar_t x, ulong y, const gr_ctx_t ctx)
{
    qqbar_mul_ui(res, x, y);
    return GR_SUCCESS;
}

int
_gr_qqbar_mul_fmpz(qqbar_t res, const qqbar_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    qqbar_mul_fmpz(res, x, y);
    return GR_SUCCESS;
}

int
_gr_qqbar_mul_fmpq(qqbar_t res, const qqbar_t x, const fmpq_t y, const gr_ctx_t ctx)
{
    qqbar_mul_fmpq(res, x, y);
    return GR_SUCCESS;
}

int
_gr_qqbar_inv(qqbar_t res, const qqbar_t x, const gr_ctx_t ctx)
{
    if (qqbar_is_zero(x))
    {
        return GR_DOMAIN;
    }
    else
    {
        qqbar_inv(res, x);
        return GR_SUCCESS;
    }
}

int
_gr_qqbar_div(qqbar_t res, const qqbar_t x, const qqbar_t y, const gr_ctx_t ctx)
{
    if (qqbar_is_zero(y))
    {
        return GR_DOMAIN;
    }
    else
    {
        qqbar_div(res, x, y);
        return GR_SUCCESS;
    }
}

int
_gr_qqbar_div_si(qqbar_t res, const qqbar_t x, slong y, const gr_ctx_t ctx)
{
    if (y == 0)
    {
        return GR_DOMAIN;
    }
    else
    {
        qqbar_div_si(res, x, y);
        return GR_SUCCESS;
    }
}

int
_gr_qqbar_div_ui(qqbar_t res, const qqbar_t x, ulong y, const gr_ctx_t ctx)
{
    if (y == 0)
    {
        return GR_DOMAIN;
    }
    else
    {
        qqbar_div_ui(res, x, y);
        return GR_SUCCESS;
    }
}

int
_gr_qqbar_div_fmpz(qqbar_t res, const qqbar_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    if (fmpz_is_zero(y))
    {
        return GR_DOMAIN;
    }
    else
    {
        qqbar_div_fmpz(res, x, y);
        return GR_SUCCESS;
    }
}

int
_gr_qqbar_div_fmpq(qqbar_t res, const qqbar_t x, const fmpq_t y, const gr_ctx_t ctx)
{
    if (fmpq_is_zero(y))
    {
        return GR_DOMAIN;
    }
    else
    {
        qqbar_div_fmpq(res, x, y);
        return GR_SUCCESS;
    }
}

truth_t
_gr_qqbar_is_invertible(const qqbar_t x, const gr_ctx_t ctx)
{
    return !qqbar_is_zero(x) ? T_TRUE : T_FALSE;
}

int
_gr_qqbar_pow_ui(qqbar_t res, const qqbar_t x, ulong exp, const gr_ctx_t ctx)
{
    qqbar_pow_ui(res, x, exp);
    return GR_SUCCESS;
}

int
_gr_qqbar_pow_si(qqbar_t res, const qqbar_t x, slong exp, const gr_ctx_t ctx)
{
    if (exp < 0 && qqbar_is_zero(x))
    {
        return GR_DOMAIN;
    }
    else
    {
        qqbar_pow_si(res, x, exp);
        return GR_SUCCESS;
    }
}

int
_gr_qqbar_pow_fmpz(qqbar_t res, const qqbar_t x, const fmpz_t exp, const gr_ctx_t ctx)
{
    if (fmpz_sgn(exp) < 0 && qqbar_is_zero(x))
    {
        return GR_DOMAIN;
    }
    else
    {
        qqbar_pow_fmpz(res, x, exp);
        return GR_SUCCESS;
    }
}

int
_gr_qqbar_pow_fmpq(qqbar_t res, const qqbar_t x, const fmpq_t exp, const gr_ctx_t ctx)
{
    if (fmpq_sgn(exp) < 0 && qqbar_is_zero(x))
    {
        return GR_DOMAIN;
    }
    else
    {
        qqbar_pow_fmpq(res, x, exp);

        /* todo: don't compute */
        if (QQBAR_CTX(ctx)->real_only && !qqbar_is_real(res))
        {
            qqbar_zero(res);
            return GR_DOMAIN;
        }
        else
        {
            return GR_SUCCESS;
        }
    }
}

int
_gr_qqbar_pow(qqbar_t res, const qqbar_t x, const qqbar_t exp, const gr_ctx_t ctx)
{
    if (qqbar_pow(res, x, exp))
    {
        if (QQBAR_CTX(ctx)->real_only && !qqbar_is_real(res))
        {
            qqbar_zero(res);
            return GR_DOMAIN;
        }
        else
        {
            return GR_SUCCESS;
        }
    }
    else
    {
        return GR_DOMAIN;
    }
}

truth_t
_gr_qqbar_is_square(const qqbar_t x, const gr_ctx_t ctx)
{
    if (QQBAR_CTX(ctx)->real_only)
        return (qqbar_sgn_re(x) >= 0) ? T_TRUE : T_FALSE;
    else
        return T_TRUE;
}

int
_gr_qqbar_sqrt(qqbar_t res, const qqbar_t x, const gr_ctx_t ctx)
{
    if (QQBAR_CTX(ctx)->real_only && qqbar_sgn_re(x) < 0)
    {
        return GR_DOMAIN;
    }
    else
    {
        qqbar_sqrt(res, x);
        return GR_SUCCESS;
    }
}

int
_gr_qqbar_rsqrt(qqbar_t res, const qqbar_t x, const gr_ctx_t ctx)
{
    if (qqbar_is_zero(x) || (QQBAR_CTX(ctx)->real_only && qqbar_sgn_re(x) < 0))
    {
        return GR_DOMAIN;
    }
    else
    {
        qqbar_rsqrt(res, x);
        return GR_SUCCESS;
    }
}

int
_gr_qqbar_numerator(qqbar_t res, const qqbar_t x, const gr_ctx_t ctx)
{
    qqbar_numerator(res, x);
    return GR_SUCCESS;
}

int
_gr_qqbar_denominator(qqbar_t res, const qqbar_t x, const gr_ctx_t ctx)
{
    fmpz_t t;
    fmpz_init(t);
    qqbar_denominator(t, x);
    qqbar_set_fmpz(res, t);
    fmpz_clear(t);
    return GR_SUCCESS;
}

/* todo: could special-case rationals */
int
_gr_qqbar_floor(qqbar_t res, const qqbar_t x, const gr_ctx_t ctx)
{
    if (qqbar_is_integer(x))
    {
        qqbar_set(res, x);
    }
    else
    {
        fmpz_t n;
        fmpz_init(n);
        qqbar_floor(n, x);
        qqbar_set_fmpz(res, n);
        fmpz_clear(n);
    }

    return GR_SUCCESS;
}

int
_gr_qqbar_ceil(qqbar_t res, const qqbar_t x, const gr_ctx_t ctx)
{
    if (qqbar_is_integer(x))
    {
        qqbar_set(res, x);
    }
    else
    {
        fmpz_t n;
        fmpz_init(n);
        qqbar_ceil(n, x);
        qqbar_set_fmpz(res, n);
        fmpz_clear(n);
    }

    return GR_SUCCESS;
}

int
_gr_qqbar_trunc(qqbar_t res, const qqbar_t x, const gr_ctx_t ctx)
{
    if (qqbar_is_integer(x))
    {
        qqbar_set(res, x);
    }
    else
    {
        int sgn = qqbar_sgn_re(x);

        if (sgn == 0)
        {
            qqbar_zero(res);
        }
        else
        {
            fmpz_t n;
            fmpz_init(n);

            if (sgn > 0)
                qqbar_floor(n, x);
            else
                qqbar_ceil(n, x);

            qqbar_set_fmpz(res, n);
            fmpz_clear(n);
        }
    }

    return GR_SUCCESS;
}

/* todo: fast numerical path */
int
_gr_qqbar_nint(qqbar_t res, const qqbar_t x, const gr_ctx_t ctx)
{
    if (qqbar_is_integer(x))
    {
        qqbar_set(res, x);
    }
    else
    {
        qqbar_t t;
        fmpz_t n;

        qqbar_init(t);
        fmpz_init(n);

        qqbar_set_d(t, 0.5);
        qqbar_add(t, x, t);
        qqbar_floor(n, t);

        if (arb_contains_int(acb_realref(QQBAR_ENCLOSURE(t))))
        {
            qqbar_re(t, t);
            if (qqbar_is_integer(t))
            {
                fmpz_t m;
                fmpz_init(m);
                qqbar_get_fmpz(m, t);
                if (fmpz_is_odd(m))
                    fmpz_sub_ui(n, n, 1);
                fmpz_clear(m);
            }
        }

        qqbar_set_fmpz(res, n);

        fmpz_clear(n);
        qqbar_clear(t);
    }

    return GR_SUCCESS;
}

int
_gr_qqbar_i(qqbar_t res, const gr_ctx_t ctx)
{
    if (ctx->which_ring == GR_CTX_REAL_ALGEBRAIC_QQBAR)
        return GR_DOMAIN;

    qqbar_i(res);
    return GR_SUCCESS;
}

int
_gr_qqbar_abs(qqbar_t res, const qqbar_t x, const gr_ctx_t ctx)
{
    qqbar_abs(res, x);
    return GR_SUCCESS;
}

/* todo: exploit when we know that the field is real */
int
_gr_qqbar_conj(qqbar_t res, const qqbar_t x, const gr_ctx_t ctx)
{
    qqbar_conj(res, x);
    return GR_SUCCESS;
}

/* todo: exploit when we know that the field is real */
int
_gr_qqbar_re(qqbar_t res, const qqbar_t x, const gr_ctx_t ctx)
{
    qqbar_re(res, x);
    return GR_SUCCESS;
}

/* todo: exploit when we know that the field is real */
int
_gr_qqbar_im(qqbar_t res, const qqbar_t x, const gr_ctx_t ctx)
{
    qqbar_im(res, x);
    return GR_SUCCESS;
}

int
_gr_qqbar_sgn(qqbar_t res, const qqbar_t x, const gr_ctx_t ctx)
{
    qqbar_sgn(res, x);
    return GR_SUCCESS;
}

int
_gr_qqbar_csgn(qqbar_t res, const qqbar_t x, const gr_ctx_t ctx)
{
    qqbar_set_si(res, qqbar_csgn(x));
    return GR_SUCCESS;
}

int
_gr_qqbar_cmp(int * res, const qqbar_t x, const qqbar_t y, const gr_ctx_t ctx)
{
    if (!qqbar_is_real(x) || !qqbar_is_real(y))
        return GR_DOMAIN;

    *res = qqbar_cmp_re(x, y);
    return GR_SUCCESS;
}

int
_gr_qqbar_cmpabs(int * res, const qqbar_t x, const qqbar_t y, const gr_ctx_t ctx)
{
    *res = qqbar_cmpabs(x, y);
    return GR_SUCCESS;
}

/* todo: 2 pi reduction for bignum numerators */

#define TRIG(fn, real_check) \
int \
_gr_qqbar_ ## fn(qqbar_t res, const qqbar_t x, const gr_ctx_t ctx) \
{ \
    if (!qqbar_is_rational(x)) \
    { \
        return GR_DOMAIN; \
    } \
    else if (COEFF_IS_MPZ(QQBAR_COEFFS(x)[0]) || COEFF_IS_MPZ(QQBAR_COEFFS(x)[1])) \
    { \
        return GR_UNABLE; \
    } \
    else \
    { \
        slong p = -QQBAR_COEFFS(x)[0], q = QQBAR_COEFFS(x)[1]; \
        if (q > QQBAR_CTX(ctx)->deg_limit) \
            return GR_UNABLE; \
        qqbar_ ## fn(res, p, q); \
        if (real_check && QQBAR_CTX(ctx)->real_only && !qqbar_is_real(res)) \
            return GR_DOMAIN; \
        return GR_SUCCESS; \
    } \
}

#define TRIG2(fn) \
int \
_gr_qqbar_ ## fn(qqbar_t res, const qqbar_t x, const gr_ctx_t ctx) \
{ \
    if (!qqbar_is_rational(x)) \
    { \
        return GR_DOMAIN; \
    } \
    else if (COEFF_IS_MPZ(QQBAR_COEFFS(x)[0]) || COEFF_IS_MPZ(QQBAR_COEFFS(x)[1])) \
    { \
        return GR_UNABLE; \
    } \
    else \
    { \
        slong p = -QQBAR_COEFFS(x)[0], q = QQBAR_COEFFS(x)[1]; \
        if (q > QQBAR_CTX(ctx)->deg_limit) \
            return GR_UNABLE; \
        return qqbar_ ## fn(res, p, q) ? GR_SUCCESS : GR_DOMAIN; \
    } \
}

#define TRIG3(fn) \
int \
_gr_qqbar_ ## fn(qqbar_t res, const qqbar_t x, const gr_ctx_t ctx) \
{ \
    fmpq_t t; \
    slong p; \
    ulong q; \
    if (!qqbar_ ## fn(&p, &q, x)) \
        return GR_DOMAIN; \
    *fmpq_numref(t) = p; \
    *fmpq_denref(t) = q; \
    qqbar_set_fmpq(res, t); \
    return GR_SUCCESS; \
}

TRIG(exp_pi_i, 1)
TRIG(sin_pi, 0)
TRIG(cos_pi, 0)
TRIG2(tan_pi)
TRIG2(cot_pi)
TRIG2(sec_pi)
TRIG2(csc_pi)

TRIG3(log_pi_i)
TRIG3(asin_pi)
TRIG3(acos_pi)
TRIG3(atan_pi)
TRIG3(acot_pi)
TRIG3(asec_pi)
TRIG3(acsc_pi)

/* todo: root of unity / is_root_of_unity */
/* todo: exploit when we know that the field is real */


/* todo: quickly skip nonreal roots over the real algebraic numbers */
int
_gr_qqbar_poly_roots_other(gr_vec_t roots, gr_vec_t mult, const gr_poly_t poly, gr_ctx_t other_ctx, int flags, gr_ctx_t ctx)
{
    if (poly->length == 0)
        return GR_DOMAIN;

    if (other_ctx->which_ring == GR_CTX_FMPZ)
    {
        gr_ctx_t ZZ;
        slong i, j, deg, deg2;
        qqbar_struct * croots;
        int status = GR_SUCCESS;

        deg = poly->length - 1;

        gr_ctx_init_fmpz(ZZ);

        gr_vec_set_length(roots, 0, ctx);
        gr_vec_set_length(mult, 0, ZZ);

        if (deg != 0)
        {
            fmpz_poly_factor_t fac;
            fmpz_poly_factor_init(fac);
            fmpz_poly_factor(fac, (const fmpz_poly_struct *) poly);

            for (i = 0; i < fac->num; i++)
            {
                deg2 = fmpz_poly_degree(fac->p + i);

                croots = _qqbar_vec_init(deg2);
                qqbar_roots_fmpz_poly(croots, fac->p + i, QQBAR_ROOTS_IRREDUCIBLE);

                for (j = 0; j < deg2; j++)
                {
                    fmpz m2 = fac->exp[i];

                    if (QQBAR_CTX(ctx)->real_only && !qqbar_is_real(croots + j))
                        continue;

                    GR_MUST_SUCCEED(gr_vec_append(roots, croots + j, ctx));
                    GR_MUST_SUCCEED(gr_vec_append(mult, &m2, ZZ));
                }

                _qqbar_vec_clear(croots, deg2);
            }

            fmpz_poly_factor_clear(fac);
        }

        /* todo: qqbar_cmp_root_order, but must sort exponents as well */

        gr_ctx_clear(ZZ);

        return status;
    }

    return GR_UNABLE;
}

truth_t
_gr_qqbar_ctx_is_algebraically_closed(gr_ctx_t ctx)
{
    return (ctx->which_ring == GR_CTX_COMPLEX_ALGEBRAIC_QQBAR) ? T_TRUE : T_FALSE;
}

truth_t
_gr_qqbar_ctx_is_ordered_ring(gr_ctx_t ctx)
{
    return (ctx->which_ring == GR_CTX_REAL_ALGEBRAIC_QQBAR) ? T_TRUE : T_FALSE;
}

int _qqbar_methods_initialized = 0;

gr_static_method_table _qqbar_methods;

gr_method_tab_input _qqbar_methods_input[] =
{
    {GR_METHOD_CTX_WRITE,       (gr_funcptr) _gr_qqbar_ctx_write},
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
                                (gr_funcptr) _gr_qqbar_ctx_is_algebraically_closed},
    {GR_METHOD_CTX_IS_ORDERED_RING,
                                (gr_funcptr) _gr_qqbar_ctx_is_ordered_ring},
    {GR_METHOD_CTX_IS_EXACT,    (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_CANONICAL,
                                (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_INIT,            (gr_funcptr) _gr_qqbar_init},

    {GR_METHOD_CLEAR,           (gr_funcptr) _gr_qqbar_clear},
    {GR_METHOD_SWAP,            (gr_funcptr) _gr_qqbar_swap},
    {GR_METHOD_SET_SHALLOW,     (gr_funcptr) _gr_qqbar_set_shallow},
    {GR_METHOD_RANDTEST,        (gr_funcptr) _gr_qqbar_randtest},
    {GR_METHOD_WRITE,           (gr_funcptr) _gr_qqbar_write},

    {GR_METHOD_ZERO,            (gr_funcptr) _gr_qqbar_zero},
    {GR_METHOD_ONE,             (gr_funcptr) _gr_qqbar_one},
    {GR_METHOD_IS_ZERO,         (gr_funcptr) _gr_qqbar_is_zero},
    {GR_METHOD_IS_ONE,          (gr_funcptr) _gr_qqbar_is_one},
    {GR_METHOD_IS_NEG_ONE,      (gr_funcptr) _gr_qqbar_is_neg_one},
    {GR_METHOD_EQUAL,           (gr_funcptr) _gr_qqbar_equal},

    {GR_METHOD_SET,             (gr_funcptr) _gr_qqbar_set},
    {GR_METHOD_SET_SI,          (gr_funcptr) _gr_qqbar_set_si},
    {GR_METHOD_SET_UI,          (gr_funcptr) _gr_qqbar_set_ui},
    {GR_METHOD_SET_FMPZ,        (gr_funcptr) _gr_qqbar_set_fmpz},
    {GR_METHOD_SET_FMPQ,        (gr_funcptr) _gr_qqbar_set_fmpq},
    {GR_METHOD_SET_D,           (gr_funcptr) _gr_qqbar_set_d},
    {GR_METHOD_SET_OTHER,       (gr_funcptr) _gr_qqbar_set_other},

    {GR_METHOD_GET_SI,          (gr_funcptr) _gr_qqbar_get_si},
    {GR_METHOD_GET_UI,          (gr_funcptr) _gr_qqbar_get_ui},
    {GR_METHOD_GET_FMPZ,        (gr_funcptr) _gr_qqbar_get_fmpz},
    {GR_METHOD_GET_FMPQ,        (gr_funcptr) _gr_qqbar_get_fmpq},
    {GR_METHOD_GET_D,           (gr_funcptr) _gr_qqbar_get_d},

    {GR_METHOD_GET_FEXPR,       (gr_funcptr) _gr_qqbar_get_fexpr},
    {GR_METHOD_GET_FEXPR_SERIALIZE,       (gr_funcptr) _gr_qqbar_get_fexpr_serialize},
    {GR_METHOD_SET_FEXPR,       (gr_funcptr) _gr_qqbar_set_fexpr},

    {GR_METHOD_NEG,             (gr_funcptr) _gr_qqbar_neg},

    {GR_METHOD_ADD,             (gr_funcptr) _gr_qqbar_add},
    {GR_METHOD_ADD_UI,          (gr_funcptr) _gr_qqbar_add_ui},
    {GR_METHOD_ADD_SI,          (gr_funcptr) _gr_qqbar_add_si},
    {GR_METHOD_ADD_FMPZ,        (gr_funcptr) _gr_qqbar_add_fmpz},
    {GR_METHOD_ADD_FMPQ,        (gr_funcptr) _gr_qqbar_add_fmpq},

    {GR_METHOD_SUB,             (gr_funcptr) _gr_qqbar_sub},
    {GR_METHOD_SUB_UI,          (gr_funcptr) _gr_qqbar_sub_ui},
    {GR_METHOD_SUB_SI,          (gr_funcptr) _gr_qqbar_sub_si},
    {GR_METHOD_SUB_FMPZ,        (gr_funcptr) _gr_qqbar_sub_fmpz},
    {GR_METHOD_SUB_FMPQ,        (gr_funcptr) _gr_qqbar_sub_fmpq},

    {GR_METHOD_MUL,             (gr_funcptr) _gr_qqbar_mul},
    {GR_METHOD_MUL_UI,          (gr_funcptr) _gr_qqbar_mul_ui},
    {GR_METHOD_MUL_SI,          (gr_funcptr) _gr_qqbar_mul_si},
    {GR_METHOD_MUL_FMPZ,        (gr_funcptr) _gr_qqbar_mul_fmpz},
    {GR_METHOD_MUL_FMPQ,        (gr_funcptr) _gr_qqbar_mul_fmpq},

    {GR_METHOD_DIV,             (gr_funcptr) _gr_qqbar_div},
    {GR_METHOD_DIV_UI,          (gr_funcptr) _gr_qqbar_div_ui},
    {GR_METHOD_DIV_SI,          (gr_funcptr) _gr_qqbar_div_si},
    {GR_METHOD_DIV_FMPZ,        (gr_funcptr) _gr_qqbar_div_fmpz},
    {GR_METHOD_DIV_FMPQ,        (gr_funcptr) _gr_qqbar_div_fmpq},

    {GR_METHOD_IS_INVERTIBLE,   (gr_funcptr) _gr_qqbar_is_invertible},
    {GR_METHOD_INV,             (gr_funcptr) _gr_qqbar_inv},

    {GR_METHOD_POW,             (gr_funcptr) _gr_qqbar_pow},
    {GR_METHOD_POW_UI,          (gr_funcptr) _gr_qqbar_pow_ui},
    {GR_METHOD_POW_SI,          (gr_funcptr) _gr_qqbar_pow_si},
    {GR_METHOD_POW_FMPZ,        (gr_funcptr) _gr_qqbar_pow_fmpz},
    {GR_METHOD_POW_FMPQ,        (gr_funcptr) _gr_qqbar_pow_fmpq},

    {GR_METHOD_IS_SQUARE,       (gr_funcptr) _gr_qqbar_is_square},
    {GR_METHOD_SQRT,            (gr_funcptr) _gr_qqbar_sqrt},
    {GR_METHOD_RSQRT,           (gr_funcptr) _gr_qqbar_rsqrt},

    {GR_METHOD_NUMERATOR,       (gr_funcptr) _gr_qqbar_numerator},
    {GR_METHOD_DENOMINATOR,     (gr_funcptr) _gr_qqbar_denominator},

    {GR_METHOD_FLOOR,           (gr_funcptr) _gr_qqbar_floor},
    {GR_METHOD_CEIL,            (gr_funcptr) _gr_qqbar_ceil},
    {GR_METHOD_TRUNC,           (gr_funcptr) _gr_qqbar_trunc},
    {GR_METHOD_NINT,            (gr_funcptr) _gr_qqbar_nint},

    {GR_METHOD_CMP,             (gr_funcptr) _gr_qqbar_cmp},
    {GR_METHOD_CMPABS,          (gr_funcptr) _gr_qqbar_cmpabs},

    {GR_METHOD_I,               (gr_funcptr) _gr_qqbar_i},
    {GR_METHOD_ABS,             (gr_funcptr) _gr_qqbar_abs},
    {GR_METHOD_CONJ,            (gr_funcptr) _gr_qqbar_conj},
    {GR_METHOD_RE,              (gr_funcptr) _gr_qqbar_re},
    {GR_METHOD_IM,              (gr_funcptr) _gr_qqbar_im},
    {GR_METHOD_SGN,             (gr_funcptr) _gr_qqbar_sgn},
    {GR_METHOD_CSGN,            (gr_funcptr) _gr_qqbar_csgn},

    {GR_METHOD_PI,              (gr_funcptr) gr_not_in_domain},

    {GR_METHOD_EXP_PI_I,        (gr_funcptr) _gr_qqbar_exp_pi_i},
    {GR_METHOD_SIN_PI,          (gr_funcptr) _gr_qqbar_sin_pi},
    {GR_METHOD_COS_PI,          (gr_funcptr) _gr_qqbar_cos_pi},
    {GR_METHOD_TAN_PI,          (gr_funcptr) _gr_qqbar_tan_pi},
    {GR_METHOD_COT_PI,          (gr_funcptr) _gr_qqbar_cot_pi},
    {GR_METHOD_SEC_PI,          (gr_funcptr) _gr_qqbar_sec_pi},
    {GR_METHOD_CSC_PI,          (gr_funcptr) _gr_qqbar_csc_pi},

    {GR_METHOD_LOG_PI_I,        (gr_funcptr) _gr_qqbar_log_pi_i},
    {GR_METHOD_ASIN_PI,          (gr_funcptr) _gr_qqbar_asin_pi},
    {GR_METHOD_ACOS_PI,          (gr_funcptr) _gr_qqbar_acos_pi},
    {GR_METHOD_ATAN_PI,          (gr_funcptr) _gr_qqbar_atan_pi},
    {GR_METHOD_ACOT_PI,          (gr_funcptr) _gr_qqbar_acot_pi},
    {GR_METHOD_ASEC_PI,          (gr_funcptr) _gr_qqbar_asec_pi},
    {GR_METHOD_ACSC_PI,          (gr_funcptr) _gr_qqbar_acsc_pi},

    {GR_METHOD_POLY_ROOTS_OTHER, (gr_funcptr) _gr_qqbar_poly_roots_other},

    {0,                         (gr_funcptr) NULL},
};

void
gr_ctx_init_real_qqbar(gr_ctx_t ctx)
{
    ctx->which_ring = GR_CTX_REAL_ALGEBRAIC_QQBAR;
    ctx->sizeof_elem = sizeof(qqbar_struct);
    ctx->size_limit = WORD_MAX;

    QQBAR_CTX(ctx)->real_only = 1;
    QQBAR_CTX(ctx)->deg_limit = WORD_MAX;
    QQBAR_CTX(ctx)->bits_limit = WORD_MAX;

    ctx->methods = _qqbar_methods;

    if (!_qqbar_methods_initialized)
    {
        gr_method_tab_init(_qqbar_methods, _qqbar_methods_input);
        _qqbar_methods_initialized = 1;
    }
}

void
gr_ctx_init_complex_qqbar(gr_ctx_t ctx)
{
    ctx->which_ring = GR_CTX_COMPLEX_ALGEBRAIC_QQBAR;
    ctx->sizeof_elem = sizeof(qqbar_struct);
    ctx->size_limit = WORD_MAX;

    QQBAR_CTX(ctx)->real_only = 0;
    QQBAR_CTX(ctx)->deg_limit = WORD_MAX;
    QQBAR_CTX(ctx)->bits_limit = WORD_MAX;

    ctx->methods = _qqbar_methods;

    if (!_qqbar_methods_initialized)
    {
        gr_method_tab_init(_qqbar_methods, _qqbar_methods_input);
        _qqbar_methods_initialized = 1;
    }
}
