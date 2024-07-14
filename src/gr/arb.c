/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"
#include "fmpz_poly_factor.h"
#include "fmpq.h"
#include "arb_poly.h"
#include "acb_poly.h"
#include "arb_mat.h"
#include "arb_fmpz_poly.h"
#include "arb_hypgeom.h"
#include "fmpzi.h"
#include "qqbar.h"
#include "gr.h"
#include "gr_generic.h"
#include "gr_vec.h"
#include "gr_poly.h"
#include "nfloat.h"

typedef struct
{
    slong prec;
}
gr_arb_ctx;

#define ARB_CTX_PREC(ring_ctx) (((gr_arb_ctx *)((ring_ctx)))->prec)

#define DEF_FUNC(fname) \
int \
_gr_arb_ ## fname(arb_t res, const arb_t x, const gr_ctx_t ctx) \
{ \
    arb_ ## fname(res, x, ARB_CTX_PREC(ctx)); \
    return GR_SUCCESS; \
} \

#define DEF_FUNC_NOPREC(fname) \
int \
_gr_arb_ ## fname(arb_t res, const arb_t x, const gr_ctx_t ctx) \
{ \
    arb_ ## fname(res, x); \
    return GR_SUCCESS; \
} \


#define DEF_2FUNC(fname) \
int \
_gr_arb_ ## fname(arb_t res1, arb_t res2, const arb_t x, const gr_ctx_t ctx) \
{ \
    arb_ ## fname(res1, res2, x, ARB_CTX_PREC(ctx)); \
    return GR_SUCCESS; \
} \

#define DEF_FUNC2(fname) \
int \
_gr_arb_ ## fname(arb_t res, const arb_t x, const arb_t y, const gr_ctx_t ctx) \
{ \
    arb_ ## fname(res, x, y, ARB_CTX_PREC(ctx)); \
    return GR_SUCCESS; \
} \

#define DEF_FUNC_SING(fname) \
int \
_gr_arb_ ## fname(arb_t res, const arb_t x, const gr_ctx_t ctx) \
{ \
    arb_ ## fname(res, x, ARB_CTX_PREC(ctx)); \
    return arb_is_finite(res) ? GR_SUCCESS : GR_UNABLE; \
} \

#define DEF_FUNC2_SING(fname) \
int \
_gr_arb_ ## fname(arb_t res, const arb_t x, const arb_t y, const gr_ctx_t ctx) \
{ \
    arb_ ## fname(res, x, y, ARB_CTX_PREC(ctx)); \
    return arb_is_finite(res) ? GR_SUCCESS : GR_UNABLE; \
} \

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

void
_gr_arb_set_shallow(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    *res = *x;
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
    /* used by polynomial printing */
    if (arb_is_exact(x))
    {
        if (arf_is_zero(arb_midref(x)))
        {
            gr_stream_write(out, "0");
            return GR_SUCCESS;
        }
        else if (arf_is_one(arb_midref(x)))
        {
            gr_stream_write(out, "1");
            return GR_SUCCESS;
        }
        else if (arf_equal_si(arb_midref(x), -1))
        {
            gr_stream_write(out, "-1");
            return GR_SUCCESS;
        }
    }

    gr_stream_write_free(out, arb_get_str(x, ARB_CTX_PREC(ctx) * 0.30102999566398 + 1, 0));
    return GR_SUCCESS;
}

int
_gr_arb_write_n(gr_stream_t out, gr_srcptr x, slong n, gr_ctx_t ctx)
{
    n = FLINT_MAX(n, 1);
    gr_stream_write_free(out, arb_get_str(x, n, ARB_STR_NO_RADIUS));
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
_gr_arb_set_str(arb_t res, const char * x, gr_ctx_t ctx)
{
    if (!arb_set_str(res, x, ARB_CTX_PREC(ctx)))
        return GR_SUCCESS;

    return gr_generic_set_str_ring_exponents(res, x, ctx);
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
_gr_arb_set_other(arb_t res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
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

        case GR_CTX_NFLOAT:
            if (NFLOAT_IS_SPECIAL(x))
            {
                if (NFLOAT_IS_ZERO(x))
                {
                    arb_zero(res);
                    return GR_SUCCESS;
                }
                else
                {
                    return GR_UNABLE;
                }
            }
            else
            {
                nfloat_get_arf(arb_midref(res), x, x_ctx);
                mag_zero(arb_radref(res));
                arb_set_round(res, res, ARB_CTX_PREC(ctx));
                return GR_SUCCESS;
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

    return gr_generic_set_other(res, x, x_ctx, ctx);
}

int
_gr_arb_set_interval_mid_rad(arb_t res, const arb_t m, const arb_t r, const gr_ctx_t ctx)
{
    mag_t rad;
    mag_init(rad);
    arb_get_mag(rad, r);
    arb_set(res, m);
    arb_add_error_mag(res, rad);
    mag_clear(rad);
    return GR_SUCCESS;
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

    if (mag_is_zero(arb_radref(x)))
        return T_FALSE;

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

DEF_FUNC_NOPREC(set)
DEF_FUNC_NOPREC(neg)
DEF_FUNC2(add)
DEF_FUNC2(sub)
DEF_FUNC2(addmul)
DEF_FUNC2(submul)
DEF_FUNC(sqr)

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
_gr_arb_mul_two(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    arb_mul_2exp_si(res, x, 1);
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

DEF_FUNC(floor)
DEF_FUNC(ceil)
DEF_FUNC(trunc)
DEF_FUNC(nint)
DEF_FUNC_NOPREC(abs)
DEF_FUNC_NOPREC(sgn)

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
_gr_arb_arg(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    if (arb_is_nonnegative(x))
    {
        arb_zero(res);
    }
    else if (arb_is_negative(x))
    {
        arb_const_pi(res, ARB_CTX_PREC(ctx));
    }
    else
    {
        arb_t t;
        arb_init(t);
        arb_const_pi(res, 2 * MAG_BITS);
        arb_union(res, res, t, ARB_CTX_PREC(ctx));
        arb_clear(t);
    }

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

DEF_FUNC(exp)
DEF_FUNC(expm1)
DEF_FUNC_SING(log1p)

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

DEF_FUNC2(min)
DEF_FUNC2(max)

DEF_FUNC(sin)
DEF_FUNC(cos)
DEF_2FUNC(sin_cos)
DEF_FUNC_SING(tan)
DEF_FUNC_SING(cot)
DEF_FUNC_SING(sec)
DEF_FUNC_SING(csc)

DEF_FUNC(sin_pi)
DEF_FUNC(cos_pi)
DEF_2FUNC(sin_cos_pi)
/* todo: detect exact singularities */
DEF_FUNC_SING(tan_pi)
DEF_FUNC_SING(cot_pi)
DEF_FUNC_SING(csc_pi)

DEF_FUNC(sinc)
DEF_FUNC(sinc_pi)

DEF_FUNC(sinh)
DEF_FUNC(cosh)
DEF_2FUNC(sinh_cosh)
DEF_FUNC(tanh)
DEF_FUNC_SING(coth)
DEF_FUNC(sech)
DEF_FUNC_SING(csch)

DEF_FUNC_SING(asin)
DEF_FUNC_SING(acos)
DEF_FUNC(atan)
DEF_FUNC2(atan2)

DEF_FUNC(asinh)
DEF_FUNC_SING(acosh)
DEF_FUNC_SING(atanh)

int
_gr_arb_lambertw(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    arb_lambertw(res, x, 0, ARB_CTX_PREC(ctx));
    return arb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int
_gr_arb_lambertw_fmpz(arb_t res, const arb_t x, const fmpz_t k, const gr_ctx_t ctx)
{
    if (fmpz_is_zero(k))
        arb_lambertw(res, x, 0, ARB_CTX_PREC(ctx));
    else if (fmpz_equal_si(k, -1))
        arb_lambertw(res, x, 1, ARB_CTX_PREC(ctx));
    else
        return GR_DOMAIN;

    return arb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
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
_gr_arb_eulernum_ui(arb_t res, ulong x, const gr_ctx_t ctx)
{
    arb_euler_number_ui(res, x, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_eulernum_fmpz(arb_t res, const fmpz_t x, const gr_ctx_t ctx)
{
    arb_euler_number_fmpz(res, x, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_bellnum_ui(arb_t res, ulong x, const gr_ctx_t ctx)
{
    arb_bell_ui(res, x, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_bellnum_fmpz(arb_t res, const fmpz_t x, const gr_ctx_t ctx)
{
    arb_bell_fmpz(res, x, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_partitions_ui(arb_t res, ulong x, const gr_ctx_t ctx)
{
    arb_partitions_ui(res, x, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_partitions_fmpz(arb_t res, const fmpz_t x, const gr_ctx_t ctx)
{
    arb_partitions_fmpz(res, x, ARB_CTX_PREC(ctx));
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
_gr_arb_erfinv(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    arb_hypgeom_erfinv(res, x, ARB_CTX_PREC(ctx));
    return arb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int
_gr_arb_erfcinv(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    arb_hypgeom_erfcinv(res, x, ARB_CTX_PREC(ctx));
    return arb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int
_gr_arb_fresnel_s(arb_t res, const arb_t x, int normalized, const gr_ctx_t ctx)
{
    arb_hypgeom_fresnel(res, NULL, x, normalized, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_fresnel_c(arb_t res, const arb_t x, int normalized, const gr_ctx_t ctx)
{
    arb_hypgeom_fresnel(NULL, res, x, normalized, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_fresnel(arb_t res1, arb_t res2, const arb_t x, int normalized, const gr_ctx_t ctx)
{
    arb_hypgeom_fresnel(res1, res2, x, normalized, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_gamma_upper(arb_t res, const arb_t x, const arb_t y, int regularized, const gr_ctx_t ctx)
{
    arb_hypgeom_gamma_upper(res, x, y, regularized, ARB_CTX_PREC(ctx));
    return arb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int
_gr_arb_gamma_lower(arb_t res, const arb_t x, const arb_t y, int regularized, const gr_ctx_t ctx)
{
    arb_hypgeom_gamma_lower(res, x, y, regularized, ARB_CTX_PREC(ctx));
    return arb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int
_gr_arb_beta_lower(arb_t res, const arb_t x, const arb_t y, const arb_t z, int regularized, const gr_ctx_t ctx)
{
    arb_hypgeom_beta_lower(res, x, y, z, regularized, ARB_CTX_PREC(ctx));
    return arb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int
_gr_arb_exp_integral(arb_t res, const arb_t x, const arb_t y, const gr_ctx_t ctx)
{
    arb_hypgeom_expint(res, x, y, ARB_CTX_PREC(ctx));
    return arb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int
_gr_arb_exp_integral_ei(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    arb_hypgeom_ei(res, x, ARB_CTX_PREC(ctx));
    return arb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int
_gr_arb_sin_integral(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    arb_hypgeom_si(res, x, ARB_CTX_PREC(ctx));
    return arb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int
_gr_arb_cos_integral(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    arb_hypgeom_ci(res, x, ARB_CTX_PREC(ctx));
    return arb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int
_gr_arb_sinh_integral(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    arb_hypgeom_shi(res, x, ARB_CTX_PREC(ctx));
    return arb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int
_gr_arb_cosh_integral(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    arb_hypgeom_chi(res, x, ARB_CTX_PREC(ctx));
    return arb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int
_gr_arb_log_integral(arb_t res, const arb_t x, int offset, const gr_ctx_t ctx)
{
    arb_hypgeom_li(res, x, offset, ARB_CTX_PREC(ctx));
    return arb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int
_gr_arb_dilog(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    arb_hypgeom_dilog(res, x, ARB_CTX_PREC(ctx));
    return arb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
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
_gr_arb_gamma_fmpz(arb_t res, const fmpz_t x, const gr_ctx_t ctx)
{
    if (fmpz_sgn(x) > 0)
    {
        arb_gamma_fmpz(res, x, ARB_CTX_PREC(ctx));
        return GR_SUCCESS;
    }
    else
    {
        return GR_DOMAIN;
    }
}

int
_gr_arb_gamma_fmpq(arb_t res, const fmpq_t x, const gr_ctx_t ctx)
{
    if (!fmpz_is_one(fmpq_denref(x)) || fmpz_sgn(fmpq_numref(x)) > 0)
    {
        arb_gamma_fmpq(res, x, ARB_CTX_PREC(ctx));
        return GR_SUCCESS;
    }
    else
    {
        return GR_DOMAIN;
    }
}

int
_gr_arb_rgamma(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    arb_rgamma(res, x, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_fac_ui(arb_t res, ulong x, const gr_ctx_t ctx)
{
    arb_fac_ui(res, x, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_fac_fmpz(arb_t res, const fmpz_t x, const gr_ctx_t ctx)
{
    int status;
    fmpz_t t;
    fmpz_init(t);
    fmpz_add_ui(t, x, 1);
    status = _gr_arb_gamma_fmpz(res, t, ctx);
    fmpz_clear(t);
    return status;
}

int
_gr_arb_rising_ui(arb_t res, const arb_t x, ulong y, const gr_ctx_t ctx)
{
    arb_rising_ui(res, x, y, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_rising(arb_t res, const arb_t x, const arb_t y, const gr_ctx_t ctx)
{
    arb_rising(res, x, y, ARB_CTX_PREC(ctx));

    if (arb_is_finite(res))
        return GR_SUCCESS;
    else
        return GR_UNABLE;
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
_gr_arb_barnes_g(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    if (arb_is_int(x) && arb_is_nonpositive(x))
    {
        return GR_DOMAIN;
    }
    else
    {
        acb_t t;
        int status;
        acb_init(t);
        acb_set_arb(t, x);
        acb_barnes_g(t, t, ARB_CTX_PREC(ctx));
        arb_swap(res, acb_realref(t));
        status = acb_is_finite(t) ? GR_SUCCESS : GR_UNABLE;
        acb_clear(t);
        return status;
    }
}

int
_gr_arb_log_barnes_g(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    if (arb_is_positive(x))
    {
        acb_t t;
        acb_init(t);
        acb_set_arb(t, x);
        acb_log_barnes_g(t, t, ARB_CTX_PREC(ctx));
        arb_swap(res, acb_realref(t));
        acb_clear(t);
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

int _gr_arb_bessel_j(arb_t res, const arb_t x, const arb_t y, const gr_ctx_t ctx) { arb_hypgeom_bessel_j(res, x, y, ARB_CTX_PREC(ctx)); return arb_is_finite(res) ? GR_SUCCESS : GR_UNABLE; }
int _gr_arb_bessel_y(arb_t res, const arb_t x, const arb_t y, const gr_ctx_t ctx) { arb_hypgeom_bessel_y(res, x, y, ARB_CTX_PREC(ctx)); return arb_is_finite(res) ? GR_SUCCESS : GR_UNABLE; }
int _gr_arb_bessel_i(arb_t res, const arb_t x, const arb_t y, const gr_ctx_t ctx) { arb_hypgeom_bessel_i(res, x, y, ARB_CTX_PREC(ctx)); return arb_is_finite(res) ? GR_SUCCESS : GR_UNABLE; }
int _gr_arb_bessel_k(arb_t res, const arb_t x, const arb_t y, const gr_ctx_t ctx) { arb_hypgeom_bessel_k(res, x, y, ARB_CTX_PREC(ctx)); return arb_is_finite(res) ? GR_SUCCESS : GR_UNABLE; }

int _gr_arb_bessel_j_y(arb_t res1, arb_t res2, const arb_t x, const arb_t y, const gr_ctx_t ctx)
{
    arb_hypgeom_bessel_jy(res1, res2, x, y, ARB_CTX_PREC(ctx));
    return (arb_is_finite(res1) && arb_is_finite(res2)) ? GR_SUCCESS : GR_UNABLE;
}

int _gr_arb_bessel_i_scaled(arb_t res, const arb_t x, const arb_t y, const gr_ctx_t ctx) { arb_hypgeom_bessel_i_scaled(res, x, y, ARB_CTX_PREC(ctx)); return arb_is_finite(res) ? GR_SUCCESS : GR_UNABLE; }
int _gr_arb_bessel_k_scaled(arb_t res, const arb_t x, const arb_t y, const gr_ctx_t ctx) { arb_hypgeom_bessel_k_scaled(res, x, y, ARB_CTX_PREC(ctx)); return arb_is_finite(res) ? GR_SUCCESS : GR_UNABLE; }

int _gr_arb_airy(arb_t res1, arb_t res2, arb_t res3, arb_t res4, const arb_t x, const gr_ctx_t ctx) { arb_hypgeom_airy(res1, res2, res3, res4, x, ARB_CTX_PREC(ctx)); return GR_SUCCESS; }
int _gr_arb_airy_ai(arb_t res, const arb_t x, const gr_ctx_t ctx) { arb_hypgeom_airy(res, NULL, NULL, NULL, x, ARB_CTX_PREC(ctx)); return GR_SUCCESS; }
int _gr_arb_airy_ai_prime(arb_t res, const arb_t x, const gr_ctx_t ctx) { arb_hypgeom_airy(NULL, res, NULL, NULL, x, ARB_CTX_PREC(ctx)); return GR_SUCCESS; }
int _gr_arb_airy_bi(arb_t res, const arb_t x, const gr_ctx_t ctx) { arb_hypgeom_airy(NULL, NULL, res, NULL, x, ARB_CTX_PREC(ctx)); return GR_SUCCESS; }
int _gr_arb_airy_bi_prime(arb_t res, const arb_t x, const gr_ctx_t ctx) { arb_hypgeom_airy(NULL, NULL, NULL, res, x, ARB_CTX_PREC(ctx)); return GR_SUCCESS; }

int _gr_arb_airy_ai_zero(arb_t res, const fmpz_t n, const gr_ctx_t ctx)
{
    if (fmpz_sgn(n) <= 0)
        return GR_DOMAIN;
    arb_hypgeom_airy_zero(res, NULL, NULL, NULL, n, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int _gr_arb_airy_ai_prime_zero(arb_t res, const fmpz_t n, const gr_ctx_t ctx)
{
    if (fmpz_sgn(n) <= 0)
        return GR_DOMAIN;
    arb_hypgeom_airy_zero(NULL, res, NULL, NULL, n, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int _gr_arb_airy_bi_zero(arb_t res, const fmpz_t n, const gr_ctx_t ctx)
{
    if (fmpz_sgn(n) <= 0)
        return GR_DOMAIN;
    arb_hypgeom_airy_zero(NULL, NULL, res, NULL, n, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int _gr_arb_airy_bi_prime_zero(arb_t res, const fmpz_t n, const gr_ctx_t ctx)
{
    if (fmpz_sgn(n) <= 0)
        return GR_DOMAIN;
    arb_hypgeom_airy_zero(NULL, NULL, NULL, res, n, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int _gr_arb_coulomb(arb_t res1, arb_t res2, arb_t res3, arb_t res4, const arb_t x, const arb_t y, const arb_t z, const gr_ctx_t ctx)
{
    /* H+, H- are typically complex */
    /* todo: document allowing NULL, or separate F+G method? */
    if (res3 == NULL && res4 == NULL)
    {
        arb_hypgeom_coulomb(res1, res2, x, y, z, ARB_CTX_PREC(ctx));
        return (arb_is_finite(res1) && arb_is_finite(res2)) ? GR_SUCCESS : GR_UNABLE;
    }
    else
    {
        return GR_UNABLE;
    }
}

int _gr_arb_coulomb_f(arb_t res, const arb_t x, const arb_t y, const arb_t z, const gr_ctx_t ctx) { arb_hypgeom_coulomb(res, NULL,  x, y, z, ARB_CTX_PREC(ctx)); return arb_is_finite(res) ? GR_SUCCESS : GR_UNABLE; }
int _gr_arb_coulomb_g(arb_t res, const arb_t x, const arb_t y, const arb_t z, const gr_ctx_t ctx) { arb_hypgeom_coulomb(NULL, res, x, y, z, ARB_CTX_PREC(ctx)); return arb_is_finite(res) ? GR_SUCCESS : GR_UNABLE; }
int _gr_arb_coulomb_hpos(arb_t res, const arb_t x, const arb_t y, const arb_t z, const gr_ctx_t ctx) { return GR_UNABLE; }
int _gr_arb_coulomb_hneg(arb_t res, const arb_t x, const arb_t y, const arb_t z, const gr_ctx_t ctx) { return GR_UNABLE; }

int _gr_arb_chebyshev_t(arb_t res, const arb_t n, const arb_t x, const gr_ctx_t ctx)
{
    arb_hypgeom_chebyshev_t(res, n, x, ARB_CTX_PREC(ctx));
    return arb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int _gr_arb_chebyshev_u(arb_t res, const arb_t n, const arb_t x, const gr_ctx_t ctx)
{
    arb_hypgeom_chebyshev_u(res, n, x, ARB_CTX_PREC(ctx));
    return arb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int _gr_arb_jacobi_p(arb_t res, const arb_t n, const arb_t a, const arb_t b, const arb_t x, const gr_ctx_t ctx)
{
    arb_hypgeom_jacobi_p(res, n, a, b, x, ARB_CTX_PREC(ctx));
    return arb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int _gr_arb_gegenbauer_c(arb_t res, const arb_t n, const arb_t m, const arb_t x, const gr_ctx_t ctx)
{
    arb_hypgeom_gegenbauer_c(res, n, m, x, ARB_CTX_PREC(ctx));
    return arb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int _gr_arb_laguerre_l(arb_t res, const arb_t n, const arb_t m, const arb_t x, const gr_ctx_t ctx)
{
    arb_hypgeom_laguerre_l(res, n, m, x, ARB_CTX_PREC(ctx));
    return arb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int _gr_arb_hermite_h(arb_t res, const arb_t n, const arb_t x, const gr_ctx_t ctx)
{
    arb_hypgeom_hermite_h(res, n, x, ARB_CTX_PREC(ctx));
    return arb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int _gr_arb_legendre_p(arb_t res, const arb_t n, const arb_t m, const arb_t x, int type, const gr_ctx_t ctx)
{
    arb_hypgeom_legendre_p(res, n, m, x, type, ARB_CTX_PREC(ctx));
    return arb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int _gr_arb_legendre_q(arb_t res, const arb_t n, const arb_t m, const arb_t x, int type, const gr_ctx_t ctx)
{
    arb_hypgeom_legendre_q(res, n, m, x, type, ARB_CTX_PREC(ctx));
    return arb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int _gr_arb_legendre_p_root_ui(arb_t res, arb_t res2, ulong n, ulong k, const gr_ctx_t ctx)
{
    if (k >= n)
        return GR_DOMAIN;

    arb_hypgeom_legendre_p_ui_root(res, res2, n, k, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int _gr_arb_hypgeom_0f1(arb_t res, const arb_t a, const arb_t x, int flags, const gr_ctx_t ctx)
{
    arb_hypgeom_0f1(res, a, x, flags, ARB_CTX_PREC(ctx));
    return arb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int _gr_arb_hypgeom_1f1(arb_t res, const arb_t a, const arb_t b, const arb_t x, int flags, const gr_ctx_t ctx)
{
    arb_hypgeom_1f1(res, a, b, x, flags, ARB_CTX_PREC(ctx));
    return arb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int _gr_arb_hypgeom_u(arb_t res, const arb_t a, const arb_t b, const arb_t x, int flags, const gr_ctx_t ctx)
{
    arb_hypgeom_u(res, a, b, x, ARB_CTX_PREC(ctx));
    return arb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int _gr_arb_hypgeom_2f1(arb_t res, const arb_t a, const arb_t b, const arb_t c, const arb_t x, int flags, const gr_ctx_t ctx)
{
    arb_hypgeom_2f1(res, a, b, c, x, flags, ARB_CTX_PREC(ctx));
    return arb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int _gr_arb_hypgeom_pfq(arb_t res, const gr_vec_t a, const gr_vec_t b, const arb_t x, int flags, const gr_ctx_t ctx)
{
    arb_hypgeom_pfq(res, a->entries, a->length, b->entries, b->length, x, flags, ARB_CTX_PREC(ctx));
    return arb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
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

DEF_FUNC2_SING(hurwitz_zeta)
DEF_FUNC2_SING(polylog)

DEF_FUNC2_SING(agm)

static void
arb_agm1(arb_t res, const arb_t x, slong prec)
{
    arb_t t;
    arb_init(t);
    arb_one(t);
    arb_agm(res, t, x, prec);
    arb_clear(t);
}

DEF_FUNC_SING(agm1)

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

int _gr_arb_poly_taylor_shift(arb_ptr res, arb_srcptr poly, slong len, const arb_t c, gr_ctx_t ctx);

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

int
_gr_arb_mat_exp(arb_mat_t res, const arb_mat_t x, gr_ctx_t ctx)
{
    if (x->r != x->c)
        return GR_DOMAIN;

    arb_mat_exp(res, x, ARB_CTX_PREC(ctx));
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
    {GR_METHOD_SET_SHALLOW,     (gr_funcptr) _gr_arb_set_shallow},
    {GR_METHOD_RANDTEST,        (gr_funcptr) _gr_arb_randtest},
    {GR_METHOD_WRITE,           (gr_funcptr) _gr_arb_write},
    {GR_METHOD_WRITE_N,         (gr_funcptr) _gr_arb_write_n},
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
    {GR_METHOD_SET_INTERVAL_MID_RAD,    (gr_funcptr) _gr_arb_set_interval_mid_rad},
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
    {GR_METHOD_ARG,             (gr_funcptr) _gr_arb_arg},
    {GR_METHOD_CMP,             (gr_funcptr) _gr_arb_cmp},
    {GR_METHOD_CMPABS,          (gr_funcptr) _gr_arb_cmpabs},
    {GR_METHOD_MIN,             (gr_funcptr) _gr_arb_min},
    {GR_METHOD_MAX,             (gr_funcptr) _gr_arb_max},
    {GR_METHOD_I,               (gr_funcptr) gr_not_in_domain},
    {GR_METHOD_PI,              (gr_funcptr) _gr_arb_pi},
    {GR_METHOD_EULER,           (gr_funcptr) _gr_arb_euler},
    {GR_METHOD_CATALAN,         (gr_funcptr) _gr_arb_catalan},
    {GR_METHOD_KHINCHIN,        (gr_funcptr) _gr_arb_khinchin},
    {GR_METHOD_GLAISHER,        (gr_funcptr) _gr_arb_glaisher},
    {GR_METHOD_EXP,             (gr_funcptr) _gr_arb_exp},
    {GR_METHOD_EXPM1,           (gr_funcptr) _gr_arb_expm1},
    {GR_METHOD_LOG,             (gr_funcptr) _gr_arb_log},
    {GR_METHOD_LOG1P,           (gr_funcptr) _gr_arb_log1p},
    {GR_METHOD_SIN,             (gr_funcptr) _gr_arb_sin},
    {GR_METHOD_COS,             (gr_funcptr) _gr_arb_cos},
    {GR_METHOD_SIN_COS,         (gr_funcptr) _gr_arb_sin_cos},
    {GR_METHOD_TAN,             (gr_funcptr) _gr_arb_tan},
    {GR_METHOD_COT,             (gr_funcptr) _gr_arb_cot},
    {GR_METHOD_SEC,             (gr_funcptr) _gr_arb_sec},
    {GR_METHOD_CSC,             (gr_funcptr) _gr_arb_csc},
    {GR_METHOD_SIN_PI,          (gr_funcptr) _gr_arb_sin_pi},
    {GR_METHOD_COS_PI,          (gr_funcptr) _gr_arb_cos_pi},
    {GR_METHOD_SIN_COS_PI,      (gr_funcptr) _gr_arb_sin_cos_pi},
    {GR_METHOD_TAN_PI,          (gr_funcptr) _gr_arb_tan_pi},
    {GR_METHOD_COT_PI,          (gr_funcptr) _gr_arb_cot_pi},
    {GR_METHOD_CSC_PI,          (gr_funcptr) _gr_arb_csc_pi},
    {GR_METHOD_SINC,            (gr_funcptr) _gr_arb_sinc},
    {GR_METHOD_SINC_PI,         (gr_funcptr) _gr_arb_sinc_pi},
    {GR_METHOD_ASIN,            (gr_funcptr) _gr_arb_asin},
    {GR_METHOD_ACOS,            (gr_funcptr) _gr_arb_acos},
    {GR_METHOD_ATAN,            (gr_funcptr) _gr_arb_atan},
    {GR_METHOD_ATAN2,           (gr_funcptr) _gr_arb_atan2},
    {GR_METHOD_SINH,            (gr_funcptr) _gr_arb_sinh},
    {GR_METHOD_COSH,            (gr_funcptr) _gr_arb_cosh},
    {GR_METHOD_SINH_COSH,       (gr_funcptr) _gr_arb_sinh_cosh},
    {GR_METHOD_TANH,            (gr_funcptr) _gr_arb_tanh},
    {GR_METHOD_COTH,            (gr_funcptr) _gr_arb_coth},
    {GR_METHOD_SECH,            (gr_funcptr) _gr_arb_sech},
    {GR_METHOD_CSCH,            (gr_funcptr) _gr_arb_csch},
    {GR_METHOD_ATANH,           (gr_funcptr) _gr_arb_atanh},
    {GR_METHOD_ASINH,           (gr_funcptr) _gr_arb_asinh},
    {GR_METHOD_ACOSH,           (gr_funcptr) _gr_arb_acosh},
    {GR_METHOD_LAMBERTW,        (gr_funcptr) _gr_arb_lambertw},
    {GR_METHOD_LAMBERTW_FMPZ,   (gr_funcptr) _gr_arb_lambertw_fmpz},
    {GR_METHOD_FAC_UI,          (gr_funcptr) _gr_arb_fac_ui},
    {GR_METHOD_FAC_FMPZ,        (gr_funcptr) _gr_arb_fac_fmpz},
    {GR_METHOD_RISING_UI,       (gr_funcptr) _gr_arb_rising_ui},
    {GR_METHOD_RISING,          (gr_funcptr) _gr_arb_rising},
    {GR_METHOD_GAMMA,           (gr_funcptr) _gr_arb_gamma},
    {GR_METHOD_GAMMA_FMPZ,      (gr_funcptr) _gr_arb_gamma_fmpz},
    {GR_METHOD_GAMMA_FMPQ,      (gr_funcptr) _gr_arb_gamma_fmpq},
    {GR_METHOD_RGAMMA,          (gr_funcptr) _gr_arb_rgamma},
    {GR_METHOD_LGAMMA,          (gr_funcptr) _gr_arb_lgamma},
    {GR_METHOD_DIGAMMA,         (gr_funcptr) _gr_arb_digamma},
    {GR_METHOD_BARNES_G,        (gr_funcptr) _gr_arb_barnes_g},
    {GR_METHOD_LOG_BARNES_G,    (gr_funcptr) _gr_arb_log_barnes_g},
    {GR_METHOD_BERNOULLI_UI,    (gr_funcptr) _gr_arb_bernoulli_ui},
    {GR_METHOD_BERNOULLI_FMPZ,  (gr_funcptr) _gr_arb_bernoulli_fmpz},
    {GR_METHOD_EULERNUM_UI,     (gr_funcptr) _gr_arb_eulernum_ui},
    {GR_METHOD_EULERNUM_FMPZ,   (gr_funcptr) _gr_arb_eulernum_fmpz},
    {GR_METHOD_BELLNUM_UI,      (gr_funcptr) _gr_arb_bellnum_ui},
    {GR_METHOD_BELLNUM_FMPZ,    (gr_funcptr) _gr_arb_bellnum_fmpz},
    {GR_METHOD_PARTITIONS_UI,   (gr_funcptr) _gr_arb_partitions_ui},
    {GR_METHOD_PARTITIONS_FMPZ, (gr_funcptr) _gr_arb_partitions_fmpz},
    {GR_METHOD_ERF,             (gr_funcptr) _gr_arb_erf},
    {GR_METHOD_ERFI,            (gr_funcptr) _gr_arb_erfi},
    {GR_METHOD_ERFC,            (gr_funcptr) _gr_arb_erfc},
    {GR_METHOD_ERFINV,          (gr_funcptr) _gr_arb_erfinv},
    {GR_METHOD_ERFCINV,         (gr_funcptr) _gr_arb_erfcinv},
    {GR_METHOD_FRESNEL_C,       (gr_funcptr) _gr_arb_fresnel_c},
    {GR_METHOD_FRESNEL_S,       (gr_funcptr) _gr_arb_fresnel_s},
    {GR_METHOD_FRESNEL,         (gr_funcptr) _gr_arb_fresnel},
    {GR_METHOD_GAMMA_UPPER,     (gr_funcptr) _gr_arb_gamma_upper},
    {GR_METHOD_GAMMA_LOWER,     (gr_funcptr) _gr_arb_gamma_lower},
    {GR_METHOD_BETA_LOWER,      (gr_funcptr) _gr_arb_beta_lower},
    {GR_METHOD_EXP_INTEGRAL,    (gr_funcptr) _gr_arb_exp_integral},
    {GR_METHOD_EXP_INTEGRAL_EI, (gr_funcptr) _gr_arb_exp_integral_ei},
    {GR_METHOD_SIN_INTEGRAL,    (gr_funcptr) _gr_arb_sin_integral},
    {GR_METHOD_COS_INTEGRAL,    (gr_funcptr) _gr_arb_cos_integral},
    {GR_METHOD_SINH_INTEGRAL,   (gr_funcptr) _gr_arb_sinh_integral},
    {GR_METHOD_COSH_INTEGRAL,   (gr_funcptr) _gr_arb_cosh_integral},
    {GR_METHOD_LOG_INTEGRAL,    (gr_funcptr) _gr_arb_log_integral},
    {GR_METHOD_DILOG,           (gr_funcptr) _gr_arb_dilog},
    {GR_METHOD_BESSEL_J,             (gr_funcptr) _gr_arb_bessel_j},
    {GR_METHOD_BESSEL_Y,             (gr_funcptr) _gr_arb_bessel_y},
    {GR_METHOD_BESSEL_I,             (gr_funcptr) _gr_arb_bessel_i},
    {GR_METHOD_BESSEL_K,             (gr_funcptr) _gr_arb_bessel_k},
    {GR_METHOD_BESSEL_J_Y,           (gr_funcptr) _gr_arb_bessel_j_y},
    {GR_METHOD_BESSEL_I_SCALED,      (gr_funcptr) _gr_arb_bessel_i_scaled},
    {GR_METHOD_BESSEL_K_SCALED,      (gr_funcptr) _gr_arb_bessel_k_scaled},
    {GR_METHOD_AIRY,                 (gr_funcptr) _gr_arb_airy},
    {GR_METHOD_AIRY_AI,              (gr_funcptr) _gr_arb_airy_ai},
    {GR_METHOD_AIRY_BI,              (gr_funcptr) _gr_arb_airy_bi},
    {GR_METHOD_AIRY_AI_PRIME,        (gr_funcptr) _gr_arb_airy_ai_prime},
    {GR_METHOD_AIRY_BI_PRIME,        (gr_funcptr) _gr_arb_airy_bi_prime},
    {GR_METHOD_AIRY_AI_ZERO,         (gr_funcptr) _gr_arb_airy_ai_zero},
    {GR_METHOD_AIRY_BI_ZERO,         (gr_funcptr) _gr_arb_airy_bi_zero},
    {GR_METHOD_AIRY_AI_PRIME_ZERO,   (gr_funcptr) _gr_arb_airy_ai_prime_zero},
    {GR_METHOD_AIRY_BI_PRIME_ZERO,   (gr_funcptr) _gr_arb_airy_bi_prime_zero},
    {GR_METHOD_COULOMB,              (gr_funcptr) _gr_arb_coulomb},
    {GR_METHOD_COULOMB_F,            (gr_funcptr) _gr_arb_coulomb_f},
    {GR_METHOD_COULOMB_G,            (gr_funcptr) _gr_arb_coulomb_g},
    {GR_METHOD_COULOMB_HNEG,         (gr_funcptr) _gr_arb_coulomb_hneg},
    {GR_METHOD_COULOMB_HPOS,         (gr_funcptr) _gr_arb_coulomb_hpos},
    {GR_METHOD_CHEBYSHEV_T,          (gr_funcptr) _gr_arb_chebyshev_t},
    {GR_METHOD_CHEBYSHEV_U,          (gr_funcptr) _gr_arb_chebyshev_u},
    {GR_METHOD_JACOBI_P,             (gr_funcptr) _gr_arb_jacobi_p},
    {GR_METHOD_GEGENBAUER_C,         (gr_funcptr) _gr_arb_gegenbauer_c},
    {GR_METHOD_LAGUERRE_L,           (gr_funcptr) _gr_arb_laguerre_l},
    {GR_METHOD_HERMITE_H,            (gr_funcptr) _gr_arb_hermite_h},
    {GR_METHOD_LEGENDRE_P,           (gr_funcptr) _gr_arb_legendre_p},
    {GR_METHOD_LEGENDRE_Q,           (gr_funcptr) _gr_arb_legendre_q},
    {GR_METHOD_LEGENDRE_P_ROOT_UI,   (gr_funcptr) _gr_arb_legendre_p_root_ui},
    {GR_METHOD_HYPGEOM_0F1,          (gr_funcptr) _gr_arb_hypgeom_0f1},
    {GR_METHOD_HYPGEOM_1F1,          (gr_funcptr) _gr_arb_hypgeom_1f1},
    {GR_METHOD_HYPGEOM_U,            (gr_funcptr) _gr_arb_hypgeom_u},
    {GR_METHOD_HYPGEOM_2F1,          (gr_funcptr) _gr_arb_hypgeom_2f1},
    {GR_METHOD_HYPGEOM_PFQ,          (gr_funcptr) _gr_arb_hypgeom_pfq},
    {GR_METHOD_ZETA,            (gr_funcptr) _gr_arb_zeta},
    {GR_METHOD_POLYLOG,         (gr_funcptr) _gr_arb_polylog},
    {GR_METHOD_HURWITZ_ZETA,    (gr_funcptr) _gr_arb_hurwitz_zeta},
    {GR_METHOD_AGM,             (gr_funcptr) _gr_arb_agm},
    {GR_METHOD_AGM1,            (gr_funcptr) _gr_arb_agm1},
    {GR_METHOD_VEC_DOT,         (gr_funcptr) _gr_arb_vec_dot},
    {GR_METHOD_VEC_DOT_REV,     (gr_funcptr) _gr_arb_vec_dot_rev},
    {GR_METHOD_POLY_MULLOW,     (gr_funcptr) _gr_arb_poly_mullow},
    {GR_METHOD_POLY_TAYLOR_SHIFT,   (gr_funcptr) _gr_arb_poly_taylor_shift},
    {GR_METHOD_POLY_ROOTS,      (gr_funcptr) _gr_arb_poly_roots},
    {GR_METHOD_POLY_ROOTS_OTHER,(gr_funcptr) _gr_arb_poly_roots_other},
    {GR_METHOD_MAT_MUL,         (gr_funcptr) _gr_arb_mat_mul},
    {GR_METHOD_MAT_DET,         (gr_funcptr) _gr_arb_mat_det},
    {GR_METHOD_MAT_EXP,         (gr_funcptr) _gr_arb_mat_exp},
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
