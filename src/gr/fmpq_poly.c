/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpq.h"
#include "fmpq_vec.h"
#include "fmpq_poly.h"
#include "gr.h"
#include "gr_poly.h"

int
_gr_fmpq_poly_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    gr_stream_write(out, "Polynomials over rationals (fmpq_poly)");
    return GR_SUCCESS;
}

void
_gr_fmpq_poly_init(fmpq_poly_t x, const gr_ctx_t ctx)
{
    fmpq_poly_init(x);
}

void
_gr_fmpq_poly_clear(fmpq_poly_t x, const gr_ctx_t ctx)
{
    fmpq_poly_clear(x);
}

void
_gr_fmpq_poly_swap(fmpq_poly_t x, fmpq_poly_t y, const gr_ctx_t ctx)
{
    fmpq_poly_t t;
    *t = *x;
    *x = *y;
    *y = *t;
}

void
_gr_fmpq_poly_set_shallow(fmpq_poly_t res, const fmpq_poly_t x, const gr_ctx_t ctx)
{
    *res = *x;
}

/* todo: limits */
int
_gr_fmpq_poly_randtest(fmpq_poly_t res, flint_rand_t state, const gr_ctx_t ctx)
{
    if (n_randint(state, 10) == 0)
        fmpq_poly_randtest(res, state, 4, 100);
    else
        fmpq_poly_randtest(res, state, 5, 10);

    return GR_SUCCESS;
}

int
_gr_fmpq_poly_write(gr_stream_t out, const fmpq_poly_t x, const gr_ctx_t ctx)
{
    gr_stream_write_free(out, fmpq_poly_get_str_pretty(x, "x"));
    return GR_SUCCESS;
}

int
_gr_fmpq_poly_zero(fmpq_poly_t x, const gr_ctx_t ctx)
{
    fmpq_poly_zero(x);
    return GR_SUCCESS;
}

int
_gr_fmpq_poly_one(fmpq_poly_t x, const gr_ctx_t ctx)
{
    fmpq_poly_one(x);
    return GR_SUCCESS;
}

int
_gr_fmpq_poly_gen(fmpq_poly_t x, const gr_ctx_t ctx)
{
    fmpq_poly_zero(x);
    fmpq_poly_set_coeff_ui(x, 1, 1);
    return GR_SUCCESS;
}

int
_gr_fmpq_poly_set_si(fmpq_poly_t res, slong v, const gr_ctx_t ctx)
{
    fmpq_poly_set_si(res, v);
    return GR_SUCCESS;
}

int
_gr_fmpq_poly_set_ui(fmpq_poly_t res, ulong v, const gr_ctx_t ctx)
{
    fmpq_poly_set_ui(res, v);
    return GR_SUCCESS;
}

int
_gr_fmpq_poly_set_fmpz(fmpq_poly_t res, const fmpz_t v, const gr_ctx_t ctx)
{
    fmpq_poly_set_fmpz(res, v);
    return GR_SUCCESS;
}

int
_gr_fmpq_poly_set_other(fmpq_poly_t res, gr_srcptr x, gr_ctx_t x_ctx, const gr_ctx_t ctx)
{
    if (x_ctx->which_ring == GR_CTX_FMPZ)
    {
        fmpq_poly_set_fmpz(res, x);
        return GR_SUCCESS;
    }

    if (x_ctx->which_ring == GR_CTX_FMPZ_POLY)
    {
        fmpq_poly_set_fmpz_poly(res, x);
        return GR_SUCCESS;
    }

    if (x_ctx->which_ring == GR_CTX_FMPQ)
    {
        fmpq_poly_set_fmpq(res, x);
        return GR_SUCCESS;
    }

    if (x_ctx->which_ring == GR_CTX_FMPQ_POLY)
    {
        fmpq_poly_set(res, x);
        return GR_SUCCESS;
    }

    if (x_ctx->which_ring == GR_CTX_GR_POLY)
    {
        if (POLYNOMIAL_ELEM_CTX(x_ctx)->which_ring == GR_CTX_FMPZ)
        {
            fmpq_poly_set_fmpz_poly(res, x);
            return GR_SUCCESS;
        }
        else
        {
            gr_ctx_t QQ;
            gr_poly_t poly;
            int status;

            gr_ctx_init_fmpq(QQ);  /* no need to free */
            gr_poly_init(poly, QQ);

            status = gr_poly_set_gr_poly_other(poly, x, POLYNOMIAL_ELEM_CTX(x_ctx), QQ);

            if (status == GR_SUCCESS)
            {
                fmpq_poly_fit_length(res, poly->length);
                _fmpq_vec_get_fmpz_vec_fmpz(res->coeffs, res->den, poly->coeffs, poly->length);
                _fmpq_poly_set_length(res, poly->length);
            }

            gr_poly_clear(poly, QQ);

            return status;
        }
    }

    if (x_ctx->which_ring == GR_CTX_GR_VEC)
    {
        gr_ctx_t QQ;
        gr_poly_t tmp, poly;
        int status;

        tmp->coeffs = ((gr_vec_struct *) x)->entries;
        tmp->length = ((gr_vec_struct *) x)->length;

        gr_ctx_init_fmpq(QQ);  /* no need to free */
        gr_poly_init(poly, QQ);

        status = gr_poly_set_gr_poly_other(poly, tmp, VECTOR_CTX(x_ctx)->base_ring, QQ);

        if (status == GR_SUCCESS)
        {
            fmpq_poly_fit_length(res, poly->length);
            _fmpq_vec_get_fmpz_vec_fmpz(res->coeffs, res->den, poly->coeffs, poly->length);
            _fmpq_poly_set_length(res, poly->length);
        }

        gr_poly_clear(poly, QQ);

        return status;
    }

    return GR_UNABLE;
}

int
_gr_fmpq_poly_get_ui(ulong * res, const fmpq_poly_t x, const gr_ctx_t ctx)
{
    if (fmpq_poly_length(x) == 0)
    {
        *res = 0;
        return GR_SUCCESS;
    }

    if (fmpq_poly_length(x) == 1 && fmpz_is_one(x->den))
    {
        if (fmpz_sgn(x->coeffs) < 0 || fmpz_cmp_ui(x->coeffs, UWORD_MAX) > 0)
        {
            return GR_DOMAIN;
        }
        else
        {
            *res = fmpz_get_ui(x->coeffs);
            return GR_SUCCESS;
        }
    }

    return GR_DOMAIN;
}

int
_gr_fmpq_poly_get_si(slong * res, const fmpq_poly_t x, const gr_ctx_t ctx)
{
    if (fmpq_poly_length(x) == 0)
    {
        *res = 0;
        return GR_SUCCESS;
    }

    if (fmpq_poly_length(x) == 1 && fmpz_is_one(x->den))
    {
        if (!fmpz_fits_si(x->coeffs))
        {
            return GR_DOMAIN;
        }
        else
        {
            *res = fmpz_get_si(x->coeffs);
            return GR_SUCCESS;
        }
    }

    return GR_DOMAIN;
}

int
_gr_fmpq_poly_get_fmpz(fmpz_t res, const fmpq_poly_t x, const gr_ctx_t ctx)
{
    if (fmpq_poly_length(x) == 0)
    {
        fmpz_zero(res);
        return GR_SUCCESS;
    }

    if (fmpq_poly_length(x) == 1 && fmpz_is_one(x->den))
    {
        fmpz_set(res, x->coeffs);
        return GR_SUCCESS;
    }

    return GR_DOMAIN;
}

int
_gr_fmpq_poly_get_fmpq(fmpq_t res, const fmpq_poly_t x, const gr_ctx_t ctx)
{
    if (fmpq_poly_length(x) == 0)
    {
        fmpq_zero(res);
        return GR_SUCCESS;
    }

    if (fmpq_poly_length(x) == 1)
    {
        fmpz_set(fmpq_numref(res), x->coeffs);
        fmpz_set(fmpq_denref(res), x->den);
        return GR_SUCCESS;
    }

    return GR_DOMAIN;
}


truth_t
_gr_fmpq_poly_is_zero(const fmpq_poly_t x, const gr_ctx_t ctx)
{
    return fmpq_poly_is_zero(x) ? T_TRUE : T_FALSE;
}

truth_t
_gr_fmpq_poly_is_one(const fmpq_poly_t x, const gr_ctx_t ctx)
{
    return fmpq_poly_is_one(x) ? T_TRUE : T_FALSE;
}

truth_t
_gr_fmpq_poly_is_neg_one(const fmpq_poly_t x, const gr_ctx_t ctx)
{
    return (x->length == 1 && x->coeffs[0] == -1 && fmpz_is_one(x->den)) ? T_TRUE : T_FALSE;
}

truth_t
_gr_fmpq_poly_equal(const fmpq_poly_t x, const fmpq_poly_t y, const gr_ctx_t ctx)
{
    return fmpq_poly_equal(x, y) ? T_TRUE : T_FALSE;
}

int
_gr_fmpq_poly_set(fmpq_poly_t res, const fmpq_poly_t x, const gr_ctx_t ctx)
{
    fmpq_poly_set(res, x);
    return GR_SUCCESS;
}

int
_gr_fmpq_poly_neg(fmpq_poly_t res, const fmpq_poly_t x, const gr_ctx_t ctx)
{
    fmpq_poly_neg(res, x);
    return GR_SUCCESS;
}

int
_gr_fmpq_poly_add(fmpq_poly_t res, const fmpq_poly_t x, const fmpq_poly_t y, const gr_ctx_t ctx)
{
    fmpq_poly_add(res, x, y);
    return GR_SUCCESS;
}

int
_gr_fmpq_poly_sub(fmpq_poly_t res, const fmpq_poly_t x, const fmpq_poly_t y, const gr_ctx_t ctx)
{
    fmpq_poly_sub(res, x, y);
    return GR_SUCCESS;
}

int
_gr_fmpq_poly_mul(fmpq_poly_t res, const fmpq_poly_t x, const fmpq_poly_t y, const gr_ctx_t ctx)
{
    fmpq_poly_mul(res, x, y);
    return GR_SUCCESS;
}

int
_gr_fmpq_poly_mul_ui(fmpq_poly_t res, const fmpq_poly_t x, ulong y, const gr_ctx_t ctx)
{
    fmpq_poly_scalar_mul_ui(res, x, y);
    return GR_SUCCESS;
}

int
_gr_fmpq_poly_mul_si(fmpq_poly_t res, const fmpq_poly_t x, slong y, const gr_ctx_t ctx)
{
    fmpq_poly_scalar_mul_si(res, x, y);
    return GR_SUCCESS;
}

int
_gr_fmpq_poly_mul_fmpz(fmpq_poly_t res, const fmpq_poly_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    fmpq_poly_scalar_mul_fmpz(res, x, y);
    return GR_SUCCESS;
}

int
_gr_fmpq_poly_mul_two(fmpq_poly_t res, const fmpq_poly_t x, const gr_ctx_t ctx)
{
    /* fmpq_poly_scalar_mul_2exp(res, x, 1); */
    fmpq_poly_scalar_mul_ui(res, x, 2);
    return GR_SUCCESS;
}

int
_gr_fmpq_poly_sqr(fmpq_poly_t res, const fmpq_poly_t x, const gr_ctx_t ctx)
{
    fmpq_poly_mul(res, x, x);
    return GR_SUCCESS;
}

int
_gr_fmpq_poly_inv(fmpq_poly_t res, const fmpq_poly_t x, const gr_ctx_t ctx)
{
    if (x->length == 1)
    {
        fmpq_poly_inv(res, x);
        return GR_SUCCESS;
    }
    else
    {
        return GR_DOMAIN;
    }
}

int
_gr_fmpq_poly_div(fmpq_poly_t res, const fmpq_poly_t x, const fmpq_poly_t y, const gr_ctx_t ctx)
{
    if (fmpq_poly_is_zero(y))
    {
        return GR_DOMAIN;
    }
    else
    {
        if (fmpq_poly_divides(res, x, y))
            return GR_SUCCESS;
        else
            return GR_DOMAIN;
    }
}

int
_gr_fmpq_poly_divexact(fmpq_poly_t res, const fmpq_poly_t x, const fmpq_poly_t y, const gr_ctx_t ctx)
{
    if (fmpq_poly_is_zero(y))
    {
        return GR_DOMAIN;
    }
    else
    {
        fmpq_poly_div(res, x, y);
        return GR_SUCCESS;
    }
}

int
_gr_fmpq_poly_divexact_ui(fmpq_poly_t res, const fmpq_poly_t x, ulong y, const gr_ctx_t ctx)
{
    if (y == 0)
    {
        return GR_DOMAIN;
    }
    else
    {
        fmpq_poly_scalar_div_ui(res, x, y);
        return GR_SUCCESS;
    }
}

int
_gr_fmpq_poly_divexact_si(fmpq_poly_t res, const fmpq_poly_t x, slong y, const gr_ctx_t ctx)
{
    if (y == 0)
    {
        return GR_DOMAIN;
    }
    else
    {
        fmpq_poly_scalar_div_si(res, x, y);
        return GR_SUCCESS;
    }
}

int
_gr_fmpq_poly_divexact_fmpz(fmpq_poly_t res, const fmpq_poly_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    if (fmpz_is_zero(y))
    {
        return GR_DOMAIN;
    }
    else
    {
        fmpq_poly_scalar_div_fmpz(res, x, y);
        return GR_SUCCESS;
    }
}


truth_t
_gr_fmpq_poly_is_invertible(const fmpq_poly_t x, const gr_ctx_t ctx)
{
    return (x->length == 1) ? T_TRUE : T_FALSE;
}

/* todo: efficient algo */
truth_t
_gr_fmpq_poly_divides(const fmpq_poly_t x, const fmpq_poly_t y, const gr_ctx_t ctx)
{
    truth_t res;
    fmpq_poly_t tmp;

    fmpq_poly_init(tmp);
    res = fmpq_poly_divides(tmp, y, x) ? T_TRUE : T_FALSE;
    fmpq_poly_clear(tmp);

    return res;
}

int
_gr_fmpq_poly_euclidean_div(fmpq_poly_t res, const fmpq_poly_t x, const fmpq_poly_t y, const gr_ctx_t ctx)
{
    if (fmpq_poly_is_zero(y))
    {
        return GR_DOMAIN;
    }
    else
    {
        fmpq_poly_div(res, x, y);
        return GR_SUCCESS;
    }
}

int
_gr_fmpq_poly_euclidean_rem(fmpq_poly_t res, const fmpq_poly_t x, const fmpq_poly_t y, const gr_ctx_t ctx)
{
    if (fmpq_poly_is_zero(y))
    {
        return GR_DOMAIN;
    }
    else
    {
        fmpq_poly_rem(res, x, y);
        return GR_SUCCESS;
    }
}

int
_gr_fmpq_poly_euclidean_divrem(fmpq_poly_t res1, fmpq_poly_t res2, const fmpq_poly_t x, const fmpq_poly_t y, const gr_ctx_t ctx)
{
    if (fmpq_poly_is_zero(y))
    {
        return GR_DOMAIN;
    }
    else
    {
        fmpq_poly_divrem(res1, res2, x, y);
        return GR_SUCCESS;
    }
}

int
_gr_fmpq_poly_gcd(fmpq_poly_t res, const fmpq_poly_t x, const fmpq_poly_t y, const gr_ctx_t ctx)
{
    fmpq_poly_gcd(res, x, y);
    return GR_SUCCESS;
}

int
_gr_fmpq_poly_lcm(fmpq_poly_t res, const fmpq_poly_t x, const fmpq_poly_t y, const gr_ctx_t ctx)
{
    fmpq_poly_lcm(res, x, y);
    return GR_SUCCESS;
}

int
_gr_fmpq_poly_pow_ui(fmpq_poly_t res, const fmpq_poly_t x, ulong exp, const gr_ctx_t ctx)
{
    fmpq_poly_pow(res, x, exp);
    return GR_SUCCESS;
}

int
_gr_fmpq_poly_pow_si(fmpq_poly_t res, const fmpq_poly_t x, slong exp, const gr_ctx_t ctx)
{
    if (exp < 0)
    {
        if (x->length == 1)
        {
            fmpq_poly_inv(res, x);
            fmpq_poly_pow(res, res, -exp);
            return GR_SUCCESS;
        }

        return GR_DOMAIN;
    }
    else
    {
        return _gr_fmpq_poly_pow_ui(res, x, exp, ctx);
    }
}

int
_gr_fmpq_poly_pow_fmpz(fmpq_poly_t res, const fmpq_poly_t x, const fmpz_t exp, gr_ctx_t ctx)
{
    if (!COEFF_IS_MPZ(*exp))
    {
        return _gr_fmpq_poly_pow_si(res, x, *exp, ctx);
    }
    else if (fmpq_poly_is_one(x) || gr_is_neg_one(x, ctx) == T_TRUE)
    {
        if (fmpq_poly_is_one(x) || fmpz_is_even(exp))
            fmpq_poly_one(res);
        else
            fmpq_poly_set_si(res, -1);
        return GR_SUCCESS;
    }
    else if (fmpq_poly_is_zero(x) && fmpz_sgn(exp) > 0)
    {
        fmpq_poly_zero(res);
        return GR_SUCCESS;
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
_gr_fmpq_poly_numerator(fmpq_poly_t res, const fmpq_poly_t x, const gr_ctx_t ctx)
{
    fmpq_poly_set(res, x);
    fmpz_one(fmpq_poly_denref(res));
    return GR_SUCCESS;
}

int
_gr_fmpq_poly_denominator(fmpq_poly_t res, const fmpq_poly_t x, const gr_ctx_t ctx)
{
    fmpq_poly_set_fmpz(res, fmpq_poly_denref(x));
    return GR_SUCCESS;
}

/*
todo: fmpq_poly_sqrt, fmpq_poly_is_square

truth_t
_gr_fmpq_poly_is_square(const fmpq_poly_t x, const gr_ctx_t ctx)
{
    truth_t res;
    fmpq_poly_t tmp;
    fmpq_poly_init(tmp);
    res = fmpq_poly_sqrt(tmp, x) ? T_TRUE : T_FALSE;
    fmpq_poly_clear(tmp);
    return res;
}

int
_gr_fmpq_poly_sqrt(fmpq_poly_t res, const fmpq_poly_t x, const gr_ctx_t ctx)
{
    if (fmpq_poly_sqrt(res, x))
    {
        return GR_SUCCESS;
    }
    else
    {
        return GR_DOMAIN;
    }
}

int
_gr_fmpq_poly_rsqrt(fmpq_poly_t res, const fmpq_poly_t x, const gr_ctx_t ctx)
{
    return GR_UNABLE;
}
*/

int _fmpq_poly_methods_initialized = 0;

gr_static_method_table _fmpq_poly_methods;

gr_method_tab_input _fmpq_poly_methods_input[] =
{
    {GR_METHOD_CTX_WRITE,       (gr_funcptr) _gr_fmpq_poly_ctx_write},
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
                                (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_EXACT,    (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_CANONICAL,
                                (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_INIT,            (gr_funcptr) _gr_fmpq_poly_init},
    {GR_METHOD_CLEAR,           (gr_funcptr) _gr_fmpq_poly_clear},
    {GR_METHOD_SWAP,            (gr_funcptr) _gr_fmpq_poly_swap},
    {GR_METHOD_SET_SHALLOW,     (gr_funcptr) _gr_fmpq_poly_set_shallow},
    {GR_METHOD_RANDTEST,        (gr_funcptr) _gr_fmpq_poly_randtest},
    {GR_METHOD_WRITE,           (gr_funcptr) _gr_fmpq_poly_write},
    {GR_METHOD_ZERO,            (gr_funcptr) _gr_fmpq_poly_zero},
    {GR_METHOD_ONE,             (gr_funcptr) _gr_fmpq_poly_one},
    {GR_METHOD_GEN,             (gr_funcptr) _gr_fmpq_poly_gen},
    {GR_METHOD_IS_ZERO,         (gr_funcptr) _gr_fmpq_poly_is_zero},
    {GR_METHOD_IS_ONE,          (gr_funcptr) _gr_fmpq_poly_is_one},
    {GR_METHOD_IS_NEG_ONE,      (gr_funcptr) _gr_fmpq_poly_is_neg_one},
    {GR_METHOD_EQUAL,           (gr_funcptr) _gr_fmpq_poly_equal},
    {GR_METHOD_SET,             (gr_funcptr) _gr_fmpq_poly_set},
    {GR_METHOD_SET_SI,          (gr_funcptr) _gr_fmpq_poly_set_si},
    {GR_METHOD_SET_UI,          (gr_funcptr) _gr_fmpq_poly_set_ui},
    {GR_METHOD_SET_FMPZ,        (gr_funcptr) _gr_fmpq_poly_set_fmpz},
    {GR_METHOD_SET_OTHER,       (gr_funcptr) _gr_fmpq_poly_set_other},
/*    {GR_METHOD_SET_STR,         (gr_funcptr) _gr_fmpq_poly_set_str}, */
    {GR_METHOD_GET_UI,          (gr_funcptr) _gr_fmpq_poly_get_ui},
    {GR_METHOD_GET_SI,          (gr_funcptr) _gr_fmpq_poly_get_si},
    {GR_METHOD_GET_FMPZ,        (gr_funcptr) _gr_fmpq_poly_get_fmpz},
    {GR_METHOD_NEG,             (gr_funcptr) _gr_fmpq_poly_neg},
    {GR_METHOD_ADD,             (gr_funcptr) _gr_fmpq_poly_add},
    {GR_METHOD_SUB,             (gr_funcptr) _gr_fmpq_poly_sub},
    {GR_METHOD_MUL,             (gr_funcptr) _gr_fmpq_poly_mul},
    {GR_METHOD_MUL_UI,          (gr_funcptr) _gr_fmpq_poly_mul_ui},
    {GR_METHOD_MUL_SI,          (gr_funcptr) _gr_fmpq_poly_mul_si},
    {GR_METHOD_MUL_FMPZ,        (gr_funcptr) _gr_fmpq_poly_mul_fmpz},
/*    {GR_METHOD_ADDMUL_UI,       (gr_funcptr) _gr_fmpq_poly_addmul_ui}, */
/*    {GR_METHOD_ADDMUL_SI,       (gr_funcptr) _gr_fmpq_poly_addmul_si}, */
/*    {GR_METHOD_ADDMUL_FMPZ,     (gr_funcptr) _gr_fmpq_poly_addmul_fmpz}, */
/*    {GR_METHOD_SUBMUL_UI,       (gr_funcptr) _gr_fmpq_poly_submul_ui}, */
/*    {GR_METHOD_SUBMUL_SI,       (gr_funcptr) _gr_fmpq_poly_submul_si}, */
/*    {GR_METHOD_SUBMUL_FMPZ,     (gr_funcptr) _gr_fmpq_poly_submul_fmpz}, */
    {GR_METHOD_MUL_TWO,         (gr_funcptr) _gr_fmpq_poly_mul_two},
    {GR_METHOD_SQR,             (gr_funcptr) _gr_fmpq_poly_sqr},
    {GR_METHOD_DIV,             (gr_funcptr) _gr_fmpq_poly_div},
    {GR_METHOD_DIVEXACT,        (gr_funcptr) _gr_fmpq_poly_divexact},
    {GR_METHOD_DIV_UI,          (gr_funcptr) _gr_fmpq_poly_divexact_ui},
    {GR_METHOD_DIV_SI,          (gr_funcptr) _gr_fmpq_poly_divexact_si},
    {GR_METHOD_DIV_FMPZ,        (gr_funcptr) _gr_fmpq_poly_divexact_fmpz},
    {GR_METHOD_DIVEXACT_UI,     (gr_funcptr) _gr_fmpq_poly_divexact_ui},
    {GR_METHOD_DIVEXACT_SI,     (gr_funcptr) _gr_fmpq_poly_divexact_si},
    {GR_METHOD_DIVEXACT_FMPZ,   (gr_funcptr) _gr_fmpq_poly_divexact_fmpz},
    {GR_METHOD_IS_INVERTIBLE,   (gr_funcptr) _gr_fmpq_poly_is_invertible},
    {GR_METHOD_INV,             (gr_funcptr) _gr_fmpq_poly_inv},
    {GR_METHOD_DIVIDES,         (gr_funcptr) _gr_fmpq_poly_divides},
    {GR_METHOD_EUCLIDEAN_DIV,   (gr_funcptr) _gr_fmpq_poly_euclidean_div},
    {GR_METHOD_EUCLIDEAN_REM,   (gr_funcptr) _gr_fmpq_poly_euclidean_rem},
    {GR_METHOD_EUCLIDEAN_DIVREM,(gr_funcptr) _gr_fmpq_poly_euclidean_divrem},
    {GR_METHOD_GCD,             (gr_funcptr) _gr_fmpq_poly_gcd},
    {GR_METHOD_LCM,             (gr_funcptr) _gr_fmpq_poly_lcm},
    {GR_METHOD_POW_UI,          (gr_funcptr) _gr_fmpq_poly_pow_ui},
    {GR_METHOD_POW_SI,          (gr_funcptr) _gr_fmpq_poly_pow_si},
    {GR_METHOD_POW_FMPZ,        (gr_funcptr) _gr_fmpq_poly_pow_fmpz},
    {GR_METHOD_NUMERATOR,       (gr_funcptr) _gr_fmpq_poly_numerator},
    {GR_METHOD_DENOMINATOR,     (gr_funcptr) _gr_fmpq_poly_denominator},
/*
    {GR_METHOD_IS_SQUARE,       (gr_funcptr) _gr_fmpq_poly_is_square},
    {GR_METHOD_SQRT,            (gr_funcptr) _gr_fmpq_poly_sqrt},
    {GR_METHOD_RSQRT,           (gr_funcptr) _gr_fmpq_poly_rsqrt},
*/
    {0,                         (gr_funcptr) NULL},
};

void
gr_ctx_init_fmpq_poly(gr_ctx_t ctx)
{
    ctx->which_ring = GR_CTX_FMPQ_POLY;
    ctx->sizeof_elem = sizeof(fmpq_poly_struct);
    ctx->size_limit = WORD_MAX;

    ctx->methods = _fmpq_poly_methods;

    if (!_fmpq_poly_methods_initialized)
    {
        gr_method_tab_init(_fmpq_poly_methods, _fmpq_poly_methods_input);
        _fmpq_poly_methods_initialized = 1;
    }
}
