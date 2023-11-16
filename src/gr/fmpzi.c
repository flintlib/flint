/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "fmpq.h"
#include "fexpr.h"
#include "qqbar.h"
#include "fmpzi.h"
#include "gr.h"

int
_gr_fmpzi_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    gr_stream_write(out, "Gaussian integer ring (fmpzi)");
    return GR_SUCCESS;
}

void
_gr_fmpzi_init(fmpzi_t x, const gr_ctx_t ctx)
{
    fmpzi_init(x);
}

void
_gr_fmpzi_clear(fmpzi_t x, const gr_ctx_t ctx)
{
    fmpzi_clear(x);
}

void
_gr_fmpzi_swap(fmpzi_t x, fmpzi_t y, const gr_ctx_t ctx)
{
    fmpzi_t t;
    *t = *x;
    *x = *y;
    *y = *t;
}

void
_gr_fmpzi_set_shallow(fmpzi_t res, const fmpzi_t x, const gr_ctx_t ctx)
{
    *res = *x;
}

/* todo: limits */
int
_gr_fmpzi_randtest(fmpzi_t res, flint_rand_t state, const gr_ctx_t ctx)
{
    switch (n_randint(state, 4))
    {
        case 0:
            fmpzi_randtest(res, state, 100);
            break;
        default:
            fmpzi_randtest(res, state, 10);
    }

    return GR_SUCCESS;
}

int
_gr_fmpzi_write(gr_stream_t out, const fmpzi_t x, const gr_ctx_t ctx)
{
    if (fmpz_is_zero(fmpzi_imagref(x)))
    {
        gr_stream_write_fmpz(out, fmpzi_realref(x));
    }
    else if (fmpz_is_zero(fmpzi_realref(x)))
    {
        if (fmpz_is_one(fmpzi_imagref(x)))
            gr_stream_write(out, "I");
        else if (fmpz_equal_si(fmpzi_imagref(x), -1))
            gr_stream_write(out, "-I");
        else
        {
            gr_stream_write_fmpz(out, fmpzi_imagref(x));
            gr_stream_write(out, "*I");
        }
    }
    else
    {
        gr_stream_write(out, "(");
        gr_stream_write_fmpz(out, fmpzi_realref(x));

        if (fmpz_is_one(fmpzi_imagref(x)))
            gr_stream_write(out, "+I)");
        else if (fmpz_equal_si(fmpzi_imagref(x), -1))
            gr_stream_write(out, "-I)");
        else
        {
            if (fmpz_sgn(fmpzi_imagref(x)) > 0)
                gr_stream_write(out, "+");
            gr_stream_write_fmpz(out, fmpzi_imagref(x));
            gr_stream_write(out, "*I)");
        }
    }

    return GR_SUCCESS;
}

int
_gr_fmpzi_zero(fmpzi_t x, const gr_ctx_t ctx)
{
    fmpzi_zero(x);
    return GR_SUCCESS;
}

int
_gr_fmpzi_one(fmpzi_t x, const gr_ctx_t ctx)
{
    fmpzi_one(x);
    return GR_SUCCESS;
}

int
_gr_fmpzi_set_si(fmpzi_t res, slong v, const gr_ctx_t ctx)
{
    fmpz_set_si(fmpzi_realref(res), v);
    fmpz_zero(fmpzi_imagref(res));
    return GR_SUCCESS;
}

int
_gr_fmpzi_set_ui(fmpzi_t res, ulong v, const gr_ctx_t ctx)
{
    fmpz_set_ui(fmpzi_realref(res), v);
    fmpz_zero(fmpzi_imagref(res));
    return GR_SUCCESS;
}

int
_gr_fmpzi_set_fmpz(fmpzi_t res, const fmpz_t v, const gr_ctx_t ctx)
{
    fmpz_set(fmpzi_realref(res), v);
    fmpz_zero(fmpzi_imagref(res));
    return GR_SUCCESS;
}

int
_gr_fmpzi_set_fmpq(fmpzi_t res, const fmpq_t v, const gr_ctx_t ctx)
{
    if (fmpz_is_one(fmpq_denref(v)))
    {
        fmpz_set(fmpzi_realref(res), fmpq_numref(v));
        fmpz_zero(fmpzi_imagref(res));
        return GR_SUCCESS;
    }
    else
    {
        return GR_DOMAIN;
    }
}

int
_gr_fmpzi_set(fmpzi_t res, const fmpzi_t x, const gr_ctx_t ctx)
{
    fmpzi_set(res, x);
    return GR_SUCCESS;
}

/* todo: move */
int
fmpzi_set_qqbar(fmpzi_t res, const qqbar_t x)
{
    if (qqbar_is_integer(x))
    {
        qqbar_get_fmpz(fmpzi_realref(res), x);
        fmpz_zero(fmpzi_imagref(res));
        return 1;
    }
    else
    {
        /* a + bi ==>  X^2 - 2aX + (a^2+b^2) */
        if (qqbar_is_algebraic_integer(x) && qqbar_degree(x) == 2 &&
            fmpz_is_even(QQBAR_COEFFS(x) + 1) && fmpz_sgn(QQBAR_COEFFS(x)) > 0)
        {
            fmpz_tdiv_q_2exp(fmpzi_realref(res), QQBAR_COEFFS(x) + 1, 1);
            fmpz_neg(fmpzi_realref(res), fmpzi_realref(res));

            fmpz_mul(fmpzi_imagref(res), fmpzi_realref(res), fmpzi_realref(res));
            fmpz_sub(fmpzi_imagref(res), QQBAR_COEFFS(x), fmpzi_imagref(res));

            /* todo: combined square root / check */
            if (fmpz_is_square(fmpzi_imagref(res)))
            {
                fmpz_sqrt(fmpzi_imagref(res), fmpzi_imagref(res));
                if (qqbar_sgn_im(x) < 0)
                    fmpz_neg(fmpzi_imagref(res), fmpzi_imagref(res));

                return 1;
            }
        }

        return 0;
    }
}

int
_gr_fmpzi_set_d(fmpzi_t res, double x, const gr_ctx_t ctx)
{
    if (x != x || x == HUGE_VAL || x == -HUGE_VAL)
        return GR_DOMAIN;

    if (x != floor(x))
        return GR_DOMAIN;

    fmpz_set_d(fmpzi_realref(res), x);
    fmpz_zero(fmpzi_imagref(res));
    return GR_SUCCESS;
}

int
_gr_fmpzi_set_other(fmpzi_t res, gr_srcptr x, gr_ctx_t x_ctx, const gr_ctx_t ctx)
{
    switch (x_ctx->which_ring)
    {
        case GR_CTX_FMPZ:
            return _gr_fmpzi_set_fmpz(res, x, ctx);

        case GR_CTX_FMPQ:
            return _gr_fmpzi_set_fmpq(res, x, ctx);

        case GR_CTX_FMPZI:
            return _gr_fmpzi_set(res, x, ctx);

        case GR_CTX_REAL_ALGEBRAIC_QQBAR:
        case GR_CTX_COMPLEX_ALGEBRAIC_QQBAR:
            return fmpzi_set_qqbar(res, x) ? GR_SUCCESS : GR_DOMAIN;
    }

    return GR_UNABLE;
}

int
_gr_fmpzi_get_fmpz(fmpz_t res, const fmpzi_t x, const gr_ctx_t ctx)
{
    if (!fmpz_is_zero(fmpzi_imagref(x)))
        return GR_DOMAIN;

    fmpz_set(res, fmpzi_realref(x));
    return GR_SUCCESS;
}

int
_gr_fmpzi_get_fmpq(fmpq_t res, const fmpzi_t x, const gr_ctx_t ctx)
{
    if (!fmpz_is_zero(fmpzi_imagref(x)))
        return GR_DOMAIN;

    fmpq_set_fmpz(res, fmpzi_realref(x));
    return GR_SUCCESS;
}

/* todo: qqbar should call an fmpzi method rather than vice versa */
void qqbar_set_fmpzi(qqbar_t res, const fmpzi_t x);

int
_gr_fmpzi_get_fexpr(fexpr_t res, const fmpzi_t x, const gr_ctx_t ctx)
{
    qqbar_t t;
    int status;
    qqbar_init(t);
    qqbar_set_fmpzi(t, x);
    status = qqbar_get_fexpr_formula(res, t, QQBAR_FORMULA_GAUSSIANS) ? GR_SUCCESS : GR_UNABLE;
    qqbar_clear(t);
    return status;
}

int
_gr_fmpzi_get_ui(ulong * res, const fmpzi_t x, const gr_ctx_t ctx)
{
    if (!fmpz_is_zero(fmpzi_imagref(x)))
        return GR_DOMAIN;

    if (fmpz_sgn(fmpzi_realref(x)) < 0 || fmpz_cmp_ui(fmpzi_realref(x), UWORD_MAX) > 0)
        return GR_DOMAIN;

    *res = fmpz_get_ui(fmpzi_realref(x));
    return GR_SUCCESS;
}

int
_gr_fmpzi_get_si(slong * res, const fmpzi_t x, const gr_ctx_t ctx)
{
    if (!fmpz_is_zero(fmpzi_imagref(x)))
        return GR_DOMAIN;

    if (!fmpz_fits_si(fmpzi_realref(x)))
        return GR_DOMAIN;

    *res = fmpz_get_si(fmpzi_realref(x));
    return GR_SUCCESS;
}

int
_gr_fmpzi_get_d(double * res, const fmpzi_t x, const gr_ctx_t ctx)
{
    if (!fmpz_is_zero(fmpzi_imagref(x)))
        return GR_DOMAIN;

    *res = fmpz_get_d(fmpzi_realref(x));
    return GR_SUCCESS;
}

truth_t
_gr_fmpzi_is_zero(const fmpzi_t x, const gr_ctx_t ctx)
{
    return fmpzi_is_zero(x) ? T_TRUE : T_FALSE;
}

truth_t
_gr_fmpzi_is_one(const fmpzi_t x, const gr_ctx_t ctx)
{
    return fmpzi_is_one(x) ? T_TRUE : T_FALSE;
}

truth_t
_gr_fmpzi_is_neg_one(const fmpzi_t x, const gr_ctx_t ctx)
{
    return (fmpz_is_zero(fmpzi_imagref(x)) && *fmpzi_realref(x) == -1) ? T_TRUE : T_FALSE;
}

truth_t
_gr_fmpzi_equal(const fmpzi_t x, const fmpzi_t y, const gr_ctx_t ctx)
{
    return fmpzi_equal(x, y) ? T_TRUE : T_FALSE;
}

int
_gr_fmpzi_neg(fmpzi_t res, const fmpzi_t x, const gr_ctx_t ctx)
{
    fmpzi_neg(res, x);
    return GR_SUCCESS;
}

int
_gr_fmpzi_add(fmpzi_t res, const fmpzi_t x, const fmpzi_t y, const gr_ctx_t ctx)
{
    fmpzi_add(res, x, y);
    return GR_SUCCESS;
}

int
_gr_fmpzi_add_si(fmpzi_t res, const fmpzi_t x, slong y, const gr_ctx_t ctx)
{
    fmpz_add_si(fmpzi_realref(res), fmpzi_realref(x), y);
    fmpz_set(fmpzi_imagref(res), fmpzi_imagref(x));
    return GR_SUCCESS;
}

int
_gr_fmpzi_add_ui(fmpzi_t res, const fmpzi_t x, ulong y, const gr_ctx_t ctx)
{
    fmpz_add_ui(fmpzi_realref(res), fmpzi_realref(x), y);
    fmpz_set(fmpzi_imagref(res), fmpzi_imagref(x));
    return GR_SUCCESS;
}

int
_gr_fmpzi_add_fmpz(fmpzi_t res, const fmpzi_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    fmpz_add(fmpzi_realref(res), fmpzi_realref(x), y);
    fmpz_set(fmpzi_imagref(res), fmpzi_imagref(x));
    return GR_SUCCESS;
}

int
_gr_fmpzi_sub(fmpzi_t res, const fmpzi_t x, const fmpzi_t y, const gr_ctx_t ctx)
{
    fmpzi_sub(res, x, y);
    return GR_SUCCESS;
}

int
_gr_fmpzi_sub_si(fmpzi_t res, const fmpzi_t x, slong y, const gr_ctx_t ctx)
{
    fmpz_sub_si(fmpzi_realref(res), fmpzi_realref(x), y);
    fmpz_set(fmpzi_imagref(res), fmpzi_imagref(x));
    return GR_SUCCESS;
}

int
_gr_fmpzi_sub_ui(fmpzi_t res, const fmpzi_t x, ulong y, const gr_ctx_t ctx)
{
    fmpz_sub_ui(fmpzi_realref(res), fmpzi_realref(x), y);
    fmpz_set(fmpzi_imagref(res), fmpzi_imagref(x));
    return GR_SUCCESS;
}

int
_gr_fmpzi_sub_fmpz(fmpzi_t res, const fmpzi_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    fmpz_sub(fmpzi_realref(res), fmpzi_realref(x), y);
    fmpz_set(fmpzi_imagref(res), fmpzi_imagref(x));
    return GR_SUCCESS;
}

int
_gr_fmpzi_mul(fmpzi_t res, const fmpzi_t x, const fmpzi_t y, const gr_ctx_t ctx)
{
    fmpzi_mul(res, x, y);
    return GR_SUCCESS;
}

int
_gr_fmpzi_mul_ui(fmpzi_t res, const fmpzi_t x, ulong y, const gr_ctx_t ctx)
{
    fmpz_mul_ui(fmpzi_realref(res), fmpzi_realref(x), y);
    fmpz_mul_ui(fmpzi_imagref(res), fmpzi_imagref(x), y);
    return GR_SUCCESS;
}

int
_gr_fmpzi_mul_si(fmpzi_t res, const fmpzi_t x, slong y, const gr_ctx_t ctx)
{
    fmpz_mul_si(fmpzi_realref(res), fmpzi_realref(x), y);
    fmpz_mul_si(fmpzi_imagref(res), fmpzi_imagref(x), y);
    return GR_SUCCESS;
}

int
_gr_fmpzi_mul_fmpz(fmpzi_t res, const fmpzi_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    fmpz_mul(fmpzi_realref(res), fmpzi_realref(x), y);
    fmpz_mul(fmpzi_imagref(res), fmpzi_imagref(x), y);
    return GR_SUCCESS;
}


/* todo: addmul, submul */

int
_gr_fmpzi_mul_two(fmpzi_t res, const fmpzi_t x, const gr_ctx_t ctx)
{
    fmpz_mul_2exp(fmpzi_realref(res), fmpzi_realref(x), 1);
    fmpz_mul_2exp(fmpzi_imagref(res), fmpzi_imagref(x), 1);
    return GR_SUCCESS;
}

int
_gr_fmpzi_sqr(fmpzi_t res, const fmpzi_t x, const gr_ctx_t ctx)
{
    fmpzi_mul(res, x, x);
    return GR_SUCCESS;
}

int
_gr_fmpzi_mul_2exp_si(fmpzi_t res, const fmpzi_t x, slong y, const gr_ctx_t ctx)
{
    if (y >= 0)
    {
        fmpz_mul_2exp(fmpzi_realref(res), fmpzi_realref(x), y);
        fmpz_mul_2exp(fmpzi_imagref(res), fmpzi_imagref(x), y);
        return GR_SUCCESS;
    }
    else if (fmpzi_is_zero(x))
    {
        fmpzi_zero(res);
        return GR_SUCCESS;
    }
    else
    {
        ulong val;
        ulong v = -(ulong) y;

        if (!fmpz_is_zero(fmpzi_realref(x)))
        {
            val = fmpz_val2(fmpzi_realref(x));
            if (val < v)
                return GR_DOMAIN;
        }

        if (!fmpz_is_zero(fmpzi_imagref(x)))
        {
            val = fmpz_val2(fmpzi_imagref(x));
            if (val < v)
                return GR_DOMAIN;
        }

        fmpz_tdiv_q_2exp(fmpzi_realref(res), fmpzi_realref(x), v);
        fmpz_tdiv_q_2exp(fmpzi_imagref(res), fmpzi_imagref(x), v);
        return GR_SUCCESS;
    }
}

int
_gr_fmpzi_mul_2exp_fmpz(fmpzi_t res, const fmpzi_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    if (fmpzi_is_zero(x))
    {
        fmpzi_zero(res);
        return GR_SUCCESS;
    }

    if (COEFF_IS_MPZ(*y))
        return GR_UNABLE;

    return _gr_fmpzi_mul_2exp_si(res, x, *y, ctx);
}

int
_gr_fmpzi_inv(fmpzi_t res, const fmpzi_t x, const gr_ctx_t ctx)
{
    if (fmpzi_is_unit(x))
    {
        if (fmpz_is_zero(fmpzi_imagref(x)))
            fmpzi_set(res, x);
        else
            fmpzi_neg(res, x);

        return GR_SUCCESS;
    }
    else
    {
        return GR_DOMAIN;
    }
}

/* todo: division optimizations (real/imag values, divexact, divisibility checking...) */

int
_gr_fmpzi_div(fmpzi_t res, const fmpzi_t x, const fmpzi_t y, const gr_ctx_t ctx)
{
    if (fmpzi_is_zero(y))
    {
        return GR_DOMAIN;
    }
    else
    {
        fmpzi_t r;
        int status;
        fmpzi_init(r);

        fmpzi_divrem_approx(res, r, x, y);
        status = fmpzi_is_zero(r) ? GR_SUCCESS : GR_DOMAIN;

        fmpzi_clear(r);
        return status;
    }
}

int
_gr_fmpzi_divexact(fmpzi_t res, const fmpzi_t x, const fmpzi_t y, const gr_ctx_t ctx)
{
    if (fmpzi_is_zero(y))
    {
        return GR_DOMAIN;
    }
    else
    {
        fmpzi_divexact(res, x, y);
        return GR_SUCCESS;
    }
}

int
_gr_fmpzi_divexact_ui(fmpzi_t res, const fmpzi_t x, ulong y, const gr_ctx_t ctx)
{
    if (y == 0)
    {
        return GR_DOMAIN;
    }
    else
    {
        fmpz_divexact_ui(fmpzi_realref(res), fmpzi_realref(x), y);
        fmpz_divexact_ui(fmpzi_imagref(res), fmpzi_imagref(x), y);
        return GR_SUCCESS;
    }
}

int
_gr_fmpzi_divexact_si(fmpzi_t res, const fmpzi_t x, slong y, const gr_ctx_t ctx)
{
    if (y == 0)
    {
        return GR_DOMAIN;
    }
    else
    {
        fmpz_divexact_si(fmpzi_realref(res), fmpzi_realref(x), y);
        fmpz_divexact_si(fmpzi_imagref(res), fmpzi_imagref(x), y);
        return GR_SUCCESS;
    }
}

int
_gr_fmpzi_divexact_fmpz(fmpzi_t res, const fmpzi_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    if (fmpz_is_zero(y))
    {
        return GR_DOMAIN;
    }
    else
    {
        fmpz_divexact(fmpzi_realref(res), fmpzi_realref(x), y);
        fmpz_divexact(fmpzi_imagref(res), fmpzi_imagref(x), y);
        return GR_SUCCESS;
    }
}

truth_t
_gr_fmpzi_is_invertible(const fmpzi_t x, const gr_ctx_t ctx)
{
    return fmpzi_is_unit(x) ? T_TRUE : T_FALSE;
}

truth_t
_gr_fmpzi_divides(const fmpzi_t x, const fmpzi_t y, const gr_ctx_t ctx)
{
    fmpzi_t q, r;
    truth_t result;

    if (fmpzi_is_zero(y))
        return T_TRUE;

    if (fmpzi_is_zero(x))
        return T_FALSE;

    fmpzi_init(q);
    fmpzi_init(r);

    fmpzi_divrem_approx(q, r, y, x);
    result = fmpzi_is_zero(r) ? T_TRUE : T_FALSE;

    fmpzi_clear(q);
    fmpzi_clear(r);
    return result;
}

int
_gr_fmpzi_euclidean_div(fmpzi_t res, const fmpzi_t x, const fmpzi_t y, const gr_ctx_t ctx)
{
    if (fmpzi_is_zero(y))
    {
        return GR_DOMAIN;
    }
    else
    {
        fmpzi_t r;
        fmpzi_init(r);
        fmpzi_divrem(res, r, x, y);
        fmpzi_clear(r);
        return GR_SUCCESS;
    }
}

int
_gr_fmpzi_euclidean_rem(fmpzi_t res, const fmpzi_t x, const fmpzi_t y, const gr_ctx_t ctx)
{
    if (fmpzi_is_zero(y))
    {
        return GR_DOMAIN;
    }
    else
    {
        fmpzi_t q;
        fmpzi_init(q);
        fmpzi_divrem(q, res, x, y);
        fmpzi_clear(q);
        return GR_SUCCESS;
    }
}

int
_gr_fmpzi_euclidean_divrem(fmpzi_t res1, fmpzi_t res2, const fmpzi_t x, const fmpzi_t y, const gr_ctx_t ctx)
{
    if (fmpzi_is_zero(y))
    {
        return GR_DOMAIN;
    }
    else
    {
        fmpzi_divrem(res1, res2, x, y);
        return GR_SUCCESS;
    }
}

int
_gr_fmpzi_gcd(fmpzi_t res, const fmpzi_t x, const fmpzi_t y, const gr_ctx_t ctx)
{
    fmpzi_gcd(res, x, y);
    return GR_SUCCESS;
}

int
_gr_fmpzi_lcm(fmpzi_t res, const fmpzi_t x, const fmpzi_t y, const gr_ctx_t ctx)
{
    if (fmpzi_is_zero(x) || fmpzi_is_zero(y))
    {
        fmpzi_zero(res);
    }
    else
    {
        fmpzi_t g;
        fmpzi_init(g);

        fmpzi_gcd(g, x, y);
        fmpzi_mul(res, x, y);

        if (!fmpzi_is_one(g))
            fmpzi_divexact(res, res, g);

        /* is this what we want? */
        fmpzi_canonicalise_unit(res, res);

        fmpzi_clear(g);
    }

    return GR_SUCCESS;
}

int
_gr_fmpzi_pow_ui(fmpzi_t res, const fmpzi_t x, ulong exp, const gr_ctx_t ctx)
{
    fmpzi_pow_ui(res, x, exp);
    return GR_SUCCESS;
}

int
_gr_fmpzi_pow_si(fmpzi_t res, const fmpzi_t x, slong exp, const gr_ctx_t ctx)
{
    if (exp < 0)
    {
        if (!fmpzi_is_unit(x))
            return GR_DOMAIN;

        exp = exp & 3;
    }

    return _gr_fmpzi_pow_ui(res, x, exp, ctx);
}

int
_gr_fmpzi_pow_fmpz(fmpzi_t res, const fmpzi_t x, const fmpz_t exp, const gr_ctx_t ctx)
{
    if (!COEFF_IS_MPZ(*exp))
    {
        return _gr_fmpzi_pow_si(res, x, *exp, ctx);
    }
    else if (fmpzi_is_unit(x))
    {
        fmpzi_pow_ui(res, x, fmpz_fdiv_ui(exp, 4));
        return GR_SUCCESS;
    }
    else if (fmpzi_is_zero(x))
    {
        if (fmpz_sgn(exp) > 0)
        {
            fmpzi_zero(res);
            return GR_SUCCESS;
        }

        return GR_DOMAIN;
    }
    else if (fmpz_sgn(exp) < 0)
    {
        return GR_DOMAIN;
    }
    else
    {
        return GR_UNABLE;
    }
}

int
_gr_fmpzi_pow(fmpzi_t res, const fmpzi_t x, const fmpzi_t exp, const gr_ctx_t ctx)
{
    if (fmpz_is_zero(fmpzi_imagref(exp)))
    {
        return _gr_fmpzi_pow_fmpz(res, x, fmpzi_realref(exp), ctx);
    }
    else if (fmpzi_is_zero(x))
    {
        if (fmpz_sgn(fmpzi_realref(exp)) > 0)
        {
            fmpzi_zero(res);
            return GR_SUCCESS;
        }
        else
        {
            return GR_DOMAIN;
        }
    }
    else if (fmpzi_is_one(x))
    {
        fmpzi_one(res);
        return GR_SUCCESS;
    }
    else
    {
        return GR_DOMAIN;
    }
}

/*
int
_gr_fmpzi_pow_fmpq(fmpzi_t res, const fmpzi_t x, const fmpq_t exp, const gr_ctx_t ctx)
{
    if (fmpz_is_one(fmpq_denref(exp)))
    {
        return _gr_fmpzi_pow_fmpz(res, x, fmpq_numref(exp), ctx);
    }
    else
    {
        qqbar_t t;
        int status;
        fmpzi_get_qqbar(t, x);

        if (qqbar_pow_fmpq(t, t, exp) && fmpzi_set_qqbar(res, t))
            status = GR_SUCCESS;
        else
            status = GR_UNABLE;

        qqbar_clear(t);
        return status;
    }
}
*/

/*
truth_t
_gr_fmpzi_is_square(const fmpzi_t x, const gr_ctx_t ctx)
{
    return fmpzi_is_square(x) ? T_TRUE : T_FALSE;
}

int
_gr_fmpzi_sqrt(fmpzi_t res, const fmpzi_t x, const gr_ctx_t ctx)
{
}

int
_gr_fmpzi_rsqrt(fmpzi_t res, const fmpzi_t x, const gr_ctx_t ctx)
{
}

int
_gr_fmpzi_abs(fmpzi_t res, const fmpzi_t x, const gr_ctx_t ctx)
{
}
*/

int
_gr_fmpzi_i(fmpzi_t res, const gr_ctx_t ctx)
{
    fmpz_zero(fmpzi_realref(res));
    fmpz_one(fmpzi_imagref(res));
    return GR_SUCCESS;
}

int
_gr_fmpzi_conj(fmpzi_t res, const fmpzi_t x, const gr_ctx_t ctx)
{
    fmpz_set(fmpzi_realref(res), fmpzi_realref(x));
    fmpz_neg(fmpzi_imagref(res), fmpzi_imagref(x));
    return GR_SUCCESS;
}

int
_gr_fmpzi_re(fmpzi_t res, const fmpzi_t x, const gr_ctx_t ctx)
{
    fmpz_set(fmpzi_realref(res), fmpzi_realref(x));
    fmpz_zero(fmpzi_imagref(res));
    return GR_SUCCESS;
}

int
_gr_fmpzi_im(fmpzi_t res, const fmpzi_t x, const gr_ctx_t ctx)
{
    fmpz_set(fmpzi_realref(res), fmpzi_imagref(x));
    fmpz_zero(fmpzi_imagref(res));
    return GR_SUCCESS;
}

/*
int
_gr_fmpzi_sgn(fmpzi_t res, const fmpzi_t x, const gr_ctx_t ctx)
{
}

int
_gr_fmpzi_cmp(int * res, const fmpzi_t x, const fmpzi_t y, const gr_ctx_t ctx)
{
}

int
_gr_fmpzi_cmpabs(int * res, const fmpzi_t x, const fmpzi_t y, const gr_ctx_t ctx)
{
}
*/


int _fmpzi_methods_initialized = 0;

gr_static_method_table _fmpzi_methods;

gr_method_tab_input _fmpzi_methods_input[] =
{
    {GR_METHOD_CTX_WRITE,       (gr_funcptr) _gr_fmpzi_ctx_write},
    {GR_METHOD_CTX_IS_RING,     (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_COMMUTATIVE_RING, (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_INTEGRAL_DOMAIN,  (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_FIELD,            (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_UNIQUE_FACTORIZATION_DOMAIN,
                                (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_FINITE,
                                (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_FINITE_CHARACTERISTIC,
                                (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_ALGEBRAICALLY_CLOSED,
                                (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_ORDERED_RING,
                                (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_EXACT,    (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_CANONICAL,
                                (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_INIT,            (gr_funcptr) _gr_fmpzi_init},
    {GR_METHOD_CLEAR,           (gr_funcptr) _gr_fmpzi_clear},
    {GR_METHOD_SWAP,            (gr_funcptr) _gr_fmpzi_swap},
    {GR_METHOD_SET_SHALLOW,     (gr_funcptr) _gr_fmpzi_set_shallow},
    {GR_METHOD_RANDTEST,        (gr_funcptr) _gr_fmpzi_randtest},
    {GR_METHOD_WRITE,           (gr_funcptr) _gr_fmpzi_write},
    {GR_METHOD_ZERO,            (gr_funcptr) _gr_fmpzi_zero},
    {GR_METHOD_ONE,             (gr_funcptr) _gr_fmpzi_one},
    {GR_METHOD_IS_ZERO,         (gr_funcptr) _gr_fmpzi_is_zero},
    {GR_METHOD_IS_ONE,          (gr_funcptr) _gr_fmpzi_is_one},
    {GR_METHOD_IS_NEG_ONE,      (gr_funcptr) _gr_fmpzi_is_neg_one},
    {GR_METHOD_EQUAL,           (gr_funcptr) _gr_fmpzi_equal},
    {GR_METHOD_SET,             (gr_funcptr) _gr_fmpzi_set},
    {GR_METHOD_SET_SI,          (gr_funcptr) _gr_fmpzi_set_si},
    {GR_METHOD_SET_UI,          (gr_funcptr) _gr_fmpzi_set_ui},
    {GR_METHOD_SET_FMPZ,        (gr_funcptr) _gr_fmpzi_set_fmpz},
    {GR_METHOD_SET_FMPQ,        (gr_funcptr) _gr_fmpzi_set_fmpq},
    {GR_METHOD_SET_D,           (gr_funcptr) _gr_fmpzi_set_d},
    {GR_METHOD_SET_OTHER,       (gr_funcptr) _gr_fmpzi_set_other},
/*    {GR_METHOD_SET_STR,         (gr_funcptr) _gr_fmpzi_set_str}, */
    {GR_METHOD_GET_FMPZ,        (gr_funcptr) _gr_fmpzi_get_fmpz},
    {GR_METHOD_GET_FMPQ,        (gr_funcptr) _gr_fmpzi_get_fmpq},
    {GR_METHOD_GET_UI,          (gr_funcptr) _gr_fmpzi_get_ui},
    {GR_METHOD_GET_SI,          (gr_funcptr) _gr_fmpzi_get_si},
    {GR_METHOD_GET_D,           (gr_funcptr) _gr_fmpzi_get_d},
    {GR_METHOD_GET_FEXPR,       (gr_funcptr) _gr_fmpzi_get_fexpr},
    {GR_METHOD_NEG,             (gr_funcptr) _gr_fmpzi_neg},
    {GR_METHOD_ADD,             (gr_funcptr) _gr_fmpzi_add},
    {GR_METHOD_ADD_UI,          (gr_funcptr) _gr_fmpzi_add_ui},
    {GR_METHOD_ADD_SI,          (gr_funcptr) _gr_fmpzi_add_si},
    {GR_METHOD_ADD_FMPZ,        (gr_funcptr) _gr_fmpzi_add_fmpz},
    {GR_METHOD_SUB,             (gr_funcptr) _gr_fmpzi_sub},
    {GR_METHOD_SUB_UI,          (gr_funcptr) _gr_fmpzi_sub_ui},
    {GR_METHOD_SUB_SI,          (gr_funcptr) _gr_fmpzi_sub_si},
    {GR_METHOD_SUB_FMPZ,        (gr_funcptr) _gr_fmpzi_sub_fmpz},
    {GR_METHOD_MUL,             (gr_funcptr) _gr_fmpzi_mul},
    {GR_METHOD_MUL_UI,          (gr_funcptr) _gr_fmpzi_mul_ui},
    {GR_METHOD_MUL_SI,          (gr_funcptr) _gr_fmpzi_mul_si},
    {GR_METHOD_MUL_FMPZ,        (gr_funcptr) _gr_fmpzi_mul_fmpz},
/*
    {GR_METHOD_ADDMUL,          (gr_funcptr) _gr_fmpzi_addmul},
    {GR_METHOD_SUBMUL,          (gr_funcptr) _gr_fmpzi_submul},
*/
    {GR_METHOD_MUL_TWO,         (gr_funcptr) _gr_fmpzi_mul_two},
    {GR_METHOD_SQR,             (gr_funcptr) _gr_fmpzi_sqr},
    {GR_METHOD_MUL_2EXP_SI,     (gr_funcptr) _gr_fmpzi_mul_2exp_si},
    {GR_METHOD_MUL_2EXP_FMPZ,   (gr_funcptr) _gr_fmpzi_mul_2exp_fmpz},
    {GR_METHOD_DIV,             (gr_funcptr) _gr_fmpzi_div},
    {GR_METHOD_DIVEXACT,        (gr_funcptr) _gr_fmpzi_divexact},
    {GR_METHOD_DIVEXACT_UI,     (gr_funcptr) _gr_fmpzi_divexact_ui},
    {GR_METHOD_DIVEXACT_SI,     (gr_funcptr) _gr_fmpzi_divexact_si},
    {GR_METHOD_DIVEXACT_FMPZ,   (gr_funcptr) _gr_fmpzi_divexact_fmpz},
    {GR_METHOD_IS_INVERTIBLE,   (gr_funcptr) _gr_fmpzi_is_invertible},
    {GR_METHOD_INV,             (gr_funcptr) _gr_fmpzi_inv},
    {GR_METHOD_DIVIDES,         (gr_funcptr) _gr_fmpzi_divides},
    {GR_METHOD_EUCLIDEAN_DIV,   (gr_funcptr) _gr_fmpzi_euclidean_div},
    {GR_METHOD_EUCLIDEAN_REM,   (gr_funcptr) _gr_fmpzi_euclidean_rem},
    {GR_METHOD_EUCLIDEAN_DIVREM,(gr_funcptr) _gr_fmpzi_euclidean_divrem},
    {GR_METHOD_GCD,             (gr_funcptr) _gr_fmpzi_gcd},
    {GR_METHOD_LCM,             (gr_funcptr) _gr_fmpzi_lcm},
    {GR_METHOD_POW_UI,          (gr_funcptr) _gr_fmpzi_pow_ui},
    {GR_METHOD_POW_SI,          (gr_funcptr) _gr_fmpzi_pow_si},
    {GR_METHOD_POW_FMPZ,        (gr_funcptr) _gr_fmpzi_pow_fmpz},
/*    {GR_METHOD_POW_FMPQ,        (gr_funcptr) _gr_fmpzi_pow_fmpq}, */
    {GR_METHOD_POW,             (gr_funcptr) _gr_fmpzi_pow},
/*
    {GR_METHOD_IS_SQUARE,       (gr_funcptr) _gr_fmpzi_is_square},
    {GR_METHOD_SQRT,            (gr_funcptr) _gr_fmpzi_sqrt},
    {GR_METHOD_RSQRT,           (gr_funcptr) _gr_fmpzi_rsqrt},
    {GR_METHOD_FLOOR,           (gr_funcptr) _gr_fmpzi_floor},
    {GR_METHOD_CEIL,            (gr_funcptr) _gr_fmpzi_ceil},
    {GR_METHOD_TRUNC,           (gr_funcptr) _gr_fmpzi_trunc},
    {GR_METHOD_NINT,            (gr_funcptr) _gr_fmpzi_nint},
*/

    {GR_METHOD_GEN,             (gr_funcptr) _gr_fmpzi_i},

    {GR_METHOD_I,               (gr_funcptr) _gr_fmpzi_i},
    {GR_METHOD_PI,              (gr_funcptr) gr_not_in_domain},
/*
    {GR_METHOD_ABS,             (gr_funcptr) _gr_fmpzi_abs},
*/
    {GR_METHOD_CONJ,            (gr_funcptr) _gr_fmpzi_conj},
    {GR_METHOD_RE,              (gr_funcptr) _gr_fmpzi_re},
    {GR_METHOD_IM,              (gr_funcptr) _gr_fmpzi_im},
/*
    {GR_METHOD_SGN,             (gr_funcptr) _gr_fmpzi_sgn},
    {GR_METHOD_CSGN,            (gr_funcptr) _gr_fmpzi_csgn},
    {GR_METHOD_CMP,             (gr_funcptr) _gr_fmpzi_cmp},
    {GR_METHOD_CMPABS,          (gr_funcptr) _gr_fmpzi_cmpabs},
    {GR_METHOD_VEC_DOT,         (gr_funcptr) _gr_fmpzi_vec_dot},
    {GR_METHOD_VEC_DOT_REV,     (gr_funcptr) _gr_fmpzi_vec_dot_rev},
    {GR_METHOD_POLY_MULLOW,     (gr_funcptr) _gr_fmpzi_poly_mullow},
    {GR_METHOD_MAT_MUL,         (gr_funcptr) _gr_fmpzi_mat_mul},
    {GR_METHOD_MAT_DET,         (gr_funcptr) _gr_fmpzi_mat_det},
*/
    {0,                         (gr_funcptr) NULL},
};

void
gr_ctx_init_fmpzi(gr_ctx_t ctx)
{
    ctx->which_ring = GR_CTX_FMPZI;
    ctx->sizeof_elem = sizeof(fmpzi_struct);
    ctx->size_limit = WORD_MAX;

    ctx->methods = _fmpzi_methods;

    if (!_fmpzi_methods_initialized)
    {
        gr_method_tab_init(_fmpzi_methods, _fmpzi_methods_input);
        _fmpzi_methods_initialized = 1;
    }
}
