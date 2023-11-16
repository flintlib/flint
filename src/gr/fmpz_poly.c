/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <string.h>
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_mat.h"
#include "fmpq.h"
#include "fmpq_poly.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_poly.h"
#include "fmpz_poly_factor.h"

#define FMPZ_POLY_CTX(ctx) POLYNOMIAL_CTX(ctx)
#define FMPZ_POLY_CTX_VAR(ctx) (FMPZ_POLY_CTX(ctx)->var)

static const char * default_var = "x";

void
_gr_fmpz_poly_ctx_clear(gr_ctx_t ctx)
{
    if (FMPZ_POLY_CTX_VAR(ctx) != default_var)
        flint_free(FMPZ_POLY_CTX_VAR(ctx));
}

int _gr_fmpz_poly_ctx_set_gen_name(gr_ctx_t ctx, const char * s)
{
    slong len;
    len = strlen(s);

    if (FMPZ_POLY_CTX_VAR(ctx) == default_var)
        FMPZ_POLY_CTX_VAR(ctx) = NULL;

    FMPZ_POLY_CTX_VAR(ctx) = flint_realloc(FMPZ_POLY_CTX_VAR(ctx), len + 1);
    memcpy(FMPZ_POLY_CTX_VAR(ctx), s, len + 1);
    return GR_SUCCESS;
}

int
_gr_fmpz_poly_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    gr_stream_write(out, "Polynomials over integers (fmpz_poly)");
    return GR_SUCCESS;
}

void
_gr_fmpz_poly_init(fmpz_poly_t x, const gr_ctx_t ctx)
{
    fmpz_poly_init(x);
}

void
_gr_fmpz_poly_clear(fmpz_poly_t x, const gr_ctx_t ctx)
{
    fmpz_poly_clear(x);
}

void
_gr_fmpz_poly_swap(fmpz_poly_t x, fmpz_poly_t y, const gr_ctx_t ctx)
{
    fmpz_poly_t t;
    *t = *x;
    *x = *y;
    *y = *t;
}

void
_gr_fmpz_poly_set_shallow(fmpz_poly_t res, const fmpz_poly_t x, const gr_ctx_t ctx)
{
    *res = *x;
}

/* todo: limits */
int
_gr_fmpz_poly_randtest(fmpz_poly_t res, flint_rand_t state, const gr_ctx_t ctx)
{
    if (n_randint(state, 10) == 0)
        fmpz_poly_randtest(res, state, 4, 100);
    else
        fmpz_poly_randtest(res, state, 5, 10);

    return GR_SUCCESS;
}

int
_gr_fmpz_poly_write(gr_stream_t out, const fmpz_poly_t x, const gr_ctx_t ctx)
{
    const char * var = FMPZ_POLY_CTX_VAR(ctx);

    if (var == NULL)
        var = "x";

    gr_stream_write_free(out, fmpz_poly_get_str_pretty(x, var));
    return GR_SUCCESS;
}

int
_gr_fmpz_poly_zero(fmpz_poly_t x, const gr_ctx_t ctx)
{
    fmpz_poly_zero(x);
    return GR_SUCCESS;
}

int
_gr_fmpz_poly_one(fmpz_poly_t x, const gr_ctx_t ctx)
{
    fmpz_poly_one(x);
    return GR_SUCCESS;
}

int
_gr_fmpz_poly_gen(fmpz_poly_t x, const gr_ctx_t ctx)
{
    fmpz_poly_zero(x);
    fmpz_poly_set_coeff_ui(x, 1, 1);
    return GR_SUCCESS;
}

int
_gr_fmpz_poly_set_si(fmpz_poly_t res, slong v, const gr_ctx_t ctx)
{
    fmpz_poly_set_si(res, v);
    return GR_SUCCESS;
}

int
_gr_fmpz_poly_set_ui(fmpz_poly_t res, ulong v, const gr_ctx_t ctx)
{
    fmpz_poly_set_ui(res, v);
    return GR_SUCCESS;
}

int
_gr_fmpz_poly_set_fmpz(fmpz_poly_t res, const fmpz_t v, const gr_ctx_t ctx)
{
    fmpz_poly_set_fmpz(res, v);
    return GR_SUCCESS;
}

int
_gr_fmpz_poly_set_other(fmpz_poly_t res, gr_srcptr x, gr_ctx_t x_ctx, const gr_ctx_t ctx)
{
    if (x_ctx->which_ring == GR_CTX_FMPZ)
    {
        fmpz_poly_set_fmpz(res, x);
        return GR_SUCCESS;
    }

    if (x_ctx->which_ring == GR_CTX_FMPZ_POLY)
    {
        fmpz_poly_set(res, x);
        return GR_SUCCESS;
    }

    if (x_ctx->which_ring == GR_CTX_FMPQ)
    {
        if (fmpz_is_one(fmpq_denref((fmpq *) x)))
        {
            fmpz_poly_set_fmpz(res, fmpq_numref((fmpq *) x));
            return GR_SUCCESS;
        }

        return GR_DOMAIN;
    }

    if (x_ctx->which_ring == GR_CTX_FMPQ_POLY)
    {
        if (fmpz_is_one(((fmpq_poly_struct *) x)->den))
        {
            fmpq_poly_get_numerator(res, x);
            return GR_SUCCESS;
        }

        return GR_DOMAIN;
    }

    if (x_ctx->which_ring == GR_CTX_GR_POLY)
    {
        if (POLYNOMIAL_ELEM_CTX(x_ctx)->which_ring == GR_CTX_FMPZ)
        {
            fmpz_poly_set(res, x);
            return GR_SUCCESS;
        }
        else
        {
            gr_ctx_t ZZ;
            gr_ctx_init_fmpz(ZZ);  /* no need to free */
            return gr_poly_set_gr_poly_other((gr_poly_struct *) res, x, POLYNOMIAL_ELEM_CTX(x_ctx), ZZ);
        }
    }

    if (x_ctx->which_ring == GR_CTX_GR_VEC)
    {
        gr_ctx_t ZZ;
        gr_poly_t tmp;
        tmp->coeffs = ((gr_vec_struct *) x)->entries;
        tmp->length = ((gr_vec_struct *) x)->length;

        gr_ctx_init_fmpz(ZZ);  /* no need to free */
        return gr_poly_set_gr_poly_other((gr_poly_struct *) res, tmp, VECTOR_CTX(x_ctx)->base_ring, ZZ);
    }

    return GR_UNABLE;
}

/*
int
_gr_fmpz_poly_set_str(fmpz_poly_t res, const char * x, const gr_ctx_t ctx)
{
    if (fmpz_poly_set_str(res, x))
        return GR_DOMAIN;

    return GR_SUCCESS;
}
*/

int
_gr_fmpz_poly_get_ui(ulong * res, const fmpz_poly_t x, const gr_ctx_t ctx)
{
    if (fmpz_poly_length(x) == 0)
    {
        *res = 0;
        return GR_SUCCESS;
    }

    if (fmpz_poly_length(x) == 1)
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
_gr_fmpz_poly_get_si(slong * res, const fmpz_poly_t x, const gr_ctx_t ctx)
{
    if (fmpz_poly_length(x) == 0)
    {
        *res = 0;
        return GR_SUCCESS;
    }

    if (fmpz_poly_length(x) == 1)
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
_gr_fmpz_poly_get_fmpz(fmpz_t res, const fmpz_poly_t x, const gr_ctx_t ctx)
{
    if (fmpz_poly_length(x) == 0)
    {
        fmpz_zero(res);
        return GR_SUCCESS;
    }

    if (fmpz_poly_length(x) == 1)
    {
        fmpz_set(res, x->coeffs);
        return GR_SUCCESS;
    }

    return GR_DOMAIN;
}

int
_gr_fmpz_poly_get_fmpq(fmpq_t res, const fmpz_poly_t x, const gr_ctx_t ctx)
{
    if (fmpz_poly_length(x) == 0)
    {
        fmpq_zero(res);
        return GR_SUCCESS;
    }

    if (fmpz_poly_length(x) == 1)
    {
        fmpq_set_fmpz(res, x->coeffs);
        return GR_SUCCESS;
    }

    return GR_DOMAIN;
}


truth_t
_gr_fmpz_poly_is_zero(const fmpz_poly_t x, const gr_ctx_t ctx)
{
    return fmpz_poly_is_zero(x) ? T_TRUE : T_FALSE;
}

truth_t
_gr_fmpz_poly_is_one(const fmpz_poly_t x, const gr_ctx_t ctx)
{
    return fmpz_poly_is_one(x) ? T_TRUE : T_FALSE;
}

truth_t
_gr_fmpz_poly_is_neg_one(const fmpz_poly_t x, const gr_ctx_t ctx)
{
    return (x->length == 1 && x->coeffs[0] == -1) ? T_TRUE : T_FALSE;
}

truth_t
_gr_fmpz_poly_equal(const fmpz_poly_t x, const fmpz_poly_t y, const gr_ctx_t ctx)
{
    return fmpz_poly_equal(x, y) ? T_TRUE : T_FALSE;
}

int
_gr_fmpz_poly_set(fmpz_poly_t res, const fmpz_poly_t x, const gr_ctx_t ctx)
{
    fmpz_poly_set(res, x);
    return GR_SUCCESS;
}

int
_gr_fmpz_poly_neg(fmpz_poly_t res, const fmpz_poly_t x, const gr_ctx_t ctx)
{
    fmpz_poly_neg(res, x);
    return GR_SUCCESS;
}

int
_gr_fmpz_poly_add(fmpz_poly_t res, const fmpz_poly_t x, const fmpz_poly_t y, const gr_ctx_t ctx)
{
    fmpz_poly_add(res, x, y);
    return GR_SUCCESS;
}

int
_gr_fmpz_poly_sub(fmpz_poly_t res, const fmpz_poly_t x, const fmpz_poly_t y, const gr_ctx_t ctx)
{
    fmpz_poly_sub(res, x, y);
    return GR_SUCCESS;
}

int
_gr_fmpz_poly_mul(fmpz_poly_t res, const fmpz_poly_t x, const fmpz_poly_t y, const gr_ctx_t ctx)
{
    fmpz_poly_mul(res, x, y);
    return GR_SUCCESS;
}

int
_gr_fmpz_poly_mul_ui(fmpz_poly_t res, const fmpz_poly_t x, ulong y, const gr_ctx_t ctx)
{
    fmpz_poly_scalar_mul_ui(res, x, y);
    return GR_SUCCESS;
}

int
_gr_fmpz_poly_mul_si(fmpz_poly_t res, const fmpz_poly_t x, slong y, const gr_ctx_t ctx)
{
    fmpz_poly_scalar_mul_si(res, x, y);
    return GR_SUCCESS;
}

int
_gr_fmpz_poly_mul_fmpz(fmpz_poly_t res, const fmpz_poly_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    fmpz_poly_scalar_mul_fmpz(res, x, y);
    return GR_SUCCESS;
}

int
_gr_fmpz_poly_addmul_ui(fmpz_poly_t res, const fmpz_poly_t x, ulong y, const gr_ctx_t ctx)
{
    fmpz_poly_scalar_addmul_ui(res, x, y);
    return GR_SUCCESS;
}

int
_gr_fmpz_poly_addmul_si(fmpz_poly_t res, const fmpz_poly_t x, slong y, const gr_ctx_t ctx)
{
    fmpz_poly_scalar_addmul_si(res, x, y);
    return GR_SUCCESS;
}

int
_gr_fmpz_poly_addmul_fmpz(fmpz_poly_t res, const fmpz_poly_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    fmpz_poly_scalar_addmul_fmpz(res, x, y);
    return GR_SUCCESS;
}

/*
todo

int
_gr_fmpz_poly_submul_ui(fmpz_poly_t res, const fmpz_poly_t x, ulong y, const gr_ctx_t ctx)
{
    fmpz_poly_scalar_submul_ui(res, x, y);
    return GR_SUCCESS;
}

int
_gr_fmpz_poly_submul_si(fmpz_poly_t res, const fmpz_poly_t x, slong y, const gr_ctx_t ctx)
{
    if (y >= 0)
        fmpz_poly_scalar_addmul_ui(res, x, y);
    else
        fmpz_poly_scalar_submul_ui(res, x, -y);
    return GR_SUCCESS;
}
*/

int
_gr_fmpz_poly_submul_fmpz(fmpz_poly_t res, const fmpz_poly_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    fmpz_poly_scalar_submul_fmpz(res, x, y);
    return GR_SUCCESS;
}

int
_gr_fmpz_poly_mul_two(fmpz_poly_t res, const fmpz_poly_t x, const gr_ctx_t ctx)
{
    fmpz_poly_scalar_mul_2exp(res, x, 1);
    return GR_SUCCESS;
}

int
_gr_fmpz_poly_sqr(fmpz_poly_t res, const fmpz_poly_t x, const gr_ctx_t ctx)
{
    fmpz_poly_sqr(res, x);
    return GR_SUCCESS;
}

int
_gr_fmpz_poly_inv(fmpz_poly_t res, const fmpz_poly_t x, const gr_ctx_t ctx)
{
    if (fmpz_poly_is_unit(x))
    {
        fmpz_poly_set(res, x);
        return GR_SUCCESS;
    }
    else
    {
        return GR_DOMAIN;
    }
}

int
_gr_fmpz_poly_div(fmpz_poly_t res, const fmpz_poly_t x, const fmpz_poly_t y, const gr_ctx_t ctx)
{
    if (fmpz_poly_is_zero(y))
    {
        return GR_DOMAIN;
    }
    else
    {
        if (fmpz_poly_divides(res, x, y))
            return GR_SUCCESS;
        else
            return GR_DOMAIN;
    }
}

int
_gr_fmpz_poly_divexact(fmpz_poly_t res, const fmpz_poly_t x, const fmpz_poly_t y, const gr_ctx_t ctx)
{
    if (fmpz_poly_is_zero(y))
    {
        return GR_DOMAIN;
    }
    else
    {
        fmpz_poly_div(res, x, y);
        return GR_SUCCESS;
    }
}

int
_gr_fmpz_poly_divexact_ui(fmpz_poly_t res, const fmpz_poly_t x, ulong y, const gr_ctx_t ctx)
{
    if (y == 0)
    {
        return GR_DOMAIN;
    }
    else
    {
        fmpz_poly_scalar_divexact_ui(res, x, y);
        return GR_SUCCESS;
    }
}

int
_gr_fmpz_poly_divexact_si(fmpz_poly_t res, const fmpz_poly_t x, slong y, const gr_ctx_t ctx)
{
    if (y == 0)
    {
        return GR_DOMAIN;
    }
    else
    {
        fmpz_poly_scalar_divexact_si(res, x, y);
        return GR_SUCCESS;
    }
}

int
_gr_fmpz_poly_divexact_fmpz(fmpz_poly_t res, const fmpz_poly_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    if (fmpz_is_zero(y))
    {
        return GR_DOMAIN;
    }
    else
    {
        fmpz_poly_scalar_divexact_fmpz(res, x, y);
        return GR_SUCCESS;
    }
}


truth_t
_gr_fmpz_poly_is_invertible(const fmpz_poly_t x, const gr_ctx_t ctx)
{
    return fmpz_poly_is_unit(x) ? T_TRUE : T_FALSE;
}

/* todo: efficient algo */
truth_t
_gr_fmpz_poly_divides(const fmpz_poly_t x, const fmpz_poly_t y, const gr_ctx_t ctx)
{
    truth_t res;
    fmpz_poly_t tmp;

    if (fmpz_poly_is_zero(y))
        return T_TRUE;

    if (fmpz_poly_is_zero(x))
        return T_FALSE;

    fmpz_poly_init(tmp);
    res = fmpz_poly_divides(tmp, y, x) ? T_TRUE : T_FALSE;
    fmpz_poly_clear(tmp);

    return res;
}

int
_gr_fmpz_poly_euclidean_div(fmpz_poly_t res, const fmpz_poly_t x, const fmpz_poly_t y, const gr_ctx_t ctx)
{
    if (fmpz_poly_is_zero(y))
    {
        return GR_DOMAIN;
    }
    else
    {
        fmpz_poly_div(res, x, y);
        return GR_SUCCESS;
    }
}

int
_gr_fmpz_poly_euclidean_rem(fmpz_poly_t res, const fmpz_poly_t x, const fmpz_poly_t y, const gr_ctx_t ctx)
{
    if (fmpz_poly_is_zero(y))
    {
        return GR_DOMAIN;
    }
    else
    {
        fmpz_poly_rem(res, x, y);
        return GR_SUCCESS;
    }
}

int
_gr_fmpz_poly_euclidean_divrem(fmpz_poly_t res1, fmpz_poly_t res2, const fmpz_poly_t x, const fmpz_poly_t y, const gr_ctx_t ctx)
{
    if (fmpz_poly_is_zero(y))
    {
        return GR_DOMAIN;
    }
    else
    {
        fmpz_poly_divrem(res1, res2, x, y);
        return GR_SUCCESS;
    }
}

int
_gr_fmpz_poly_gcd(fmpz_poly_t res, const fmpz_poly_t x, const fmpz_poly_t y, const gr_ctx_t ctx)
{
    fmpz_poly_gcd(res, x, y);
    return GR_SUCCESS;
}

int
_gr_fmpz_poly_lcm(fmpz_poly_t res, const fmpz_poly_t x, const fmpz_poly_t y, const gr_ctx_t ctx)
{
    fmpz_poly_lcm(res, x, y);
    return GR_SUCCESS;
}

int
_gr_fmpz_poly_pow_ui(fmpz_poly_t res, const fmpz_poly_t x, ulong exp, const gr_ctx_t ctx)
{
    fmpz_poly_pow(res, x, exp);
    return GR_SUCCESS;
}

int
_gr_fmpz_poly_pow_si(fmpz_poly_t res, const fmpz_poly_t x, slong exp, const gr_ctx_t ctx)
{
    if (exp < 0)
    {
        if (fmpz_poly_is_unit(x))
        {
            if (fmpz_poly_is_one(x) || exp % 2 == 0)
                fmpz_poly_one(res);
            else
                fmpz_poly_set_si(res, -1);

            return GR_SUCCESS;
        }

        return GR_DOMAIN;
    }
    else
    {
        return _gr_fmpz_poly_pow_ui(res, x, exp, ctx);
    }
}

int
_gr_fmpz_poly_pow_fmpz(fmpz_poly_t res, const fmpz_poly_t x, const fmpz_t exp, const gr_ctx_t ctx)
{
    if (!COEFF_IS_MPZ(*exp))
    {
        return _gr_fmpz_poly_pow_si(res, x, *exp, ctx);
    }
    else if (fmpz_poly_is_unit(x))
    {
        if (fmpz_poly_is_one(x) || fmpz_is_even(exp))
            fmpz_poly_one(res);
        else
            fmpz_poly_set_si(res, -1);
        return GR_SUCCESS;
    }
    else if (fmpz_poly_is_zero(x) && fmpz_sgn(exp) > 0)
    {
        fmpz_poly_zero(res);
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

truth_t
_gr_fmpz_poly_is_square(const fmpz_poly_t x, const gr_ctx_t ctx)
{
    /* todo: fmpz_poly_is_square */
    truth_t res;
    fmpz_poly_t tmp;
    fmpz_poly_init(tmp);
    res = fmpz_poly_sqrt(tmp, x) ? T_TRUE : T_FALSE;
    fmpz_poly_clear(tmp);
    return res;
}

int
_gr_fmpz_poly_sqrt(fmpz_poly_t res, const fmpz_poly_t x, const gr_ctx_t ctx)
{
    if (fmpz_poly_sqrt(res, x))
    {
        return GR_SUCCESS;
    }
    else
    {
        return GR_DOMAIN;
    }
}

int
_gr_fmpz_poly_rsqrt(fmpz_poly_t res, const fmpz_poly_t x, const gr_ctx_t ctx)
{
    if (fmpz_poly_is_one(x))
    {
        fmpz_poly_one(res);
        return GR_SUCCESS;
    }
    else
    {
        return GR_DOMAIN;
    }
}

int
_gr_fmpz_poly_factor(fmpz_poly_t c, gr_vec_t factors, gr_vec_t exponents, gr_srcptr x, int flags, gr_ctx_t ctx)
{
    fmpz_poly_factor_t fac;
    gr_ctx_t ZZ;
    slong i;

    fmpz_poly_factor_init(fac);
    fmpz_poly_factor(fac, x);

    fmpz_poly_set_fmpz(c, &fac->c);

    gr_ctx_init_fmpz(ZZ);

    gr_vec_set_length(factors, fac->num, ctx);
    gr_vec_set_length(exponents, fac->num, ZZ);

    for (i = 0; i < fac->num; i++)
    {
        fmpz_poly_swap((fmpz_poly_struct *) (factors->entries) + i, fac->p + i);
        fmpz_set_ui((fmpz *) (exponents->entries) + i, fac->exp[i]);
    }

    gr_ctx_clear(ZZ);

    fmpz_poly_factor_clear(fac);

    return GR_SUCCESS;
}

int _fmpz_poly_methods_initialized = 0;

gr_static_method_table _fmpz_poly_methods;

gr_method_tab_input _fmpz_poly_methods_input[] =
{
    {GR_METHOD_CTX_WRITE,       (gr_funcptr) _gr_fmpz_poly_ctx_write},
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

    {GR_METHOD_CTX_SET_GEN_NAME,    (gr_funcptr) _gr_fmpz_poly_ctx_set_gen_name},

    {GR_METHOD_INIT,            (gr_funcptr) _gr_fmpz_poly_init},
    {GR_METHOD_CLEAR,           (gr_funcptr) _gr_fmpz_poly_clear},
    {GR_METHOD_SWAP,            (gr_funcptr) _gr_fmpz_poly_swap},
    {GR_METHOD_SET_SHALLOW,     (gr_funcptr) _gr_fmpz_poly_set_shallow},
    {GR_METHOD_RANDTEST,        (gr_funcptr) _gr_fmpz_poly_randtest},
    {GR_METHOD_WRITE,           (gr_funcptr) _gr_fmpz_poly_write},
    {GR_METHOD_ZERO,            (gr_funcptr) _gr_fmpz_poly_zero},
    {GR_METHOD_ONE,             (gr_funcptr) _gr_fmpz_poly_one},
    {GR_METHOD_GEN,             (gr_funcptr) _gr_fmpz_poly_gen},
    {GR_METHOD_IS_ZERO,         (gr_funcptr) _gr_fmpz_poly_is_zero},
    {GR_METHOD_IS_ONE,          (gr_funcptr) _gr_fmpz_poly_is_one},
    {GR_METHOD_IS_NEG_ONE,      (gr_funcptr) _gr_fmpz_poly_is_neg_one},
    {GR_METHOD_EQUAL,           (gr_funcptr) _gr_fmpz_poly_equal},
    {GR_METHOD_SET,             (gr_funcptr) _gr_fmpz_poly_set},
    {GR_METHOD_SET_SI,          (gr_funcptr) _gr_fmpz_poly_set_si},
    {GR_METHOD_SET_UI,          (gr_funcptr) _gr_fmpz_poly_set_ui},
    {GR_METHOD_SET_FMPZ,        (gr_funcptr) _gr_fmpz_poly_set_fmpz},
    {GR_METHOD_SET_OTHER,       (gr_funcptr) _gr_fmpz_poly_set_other},
/*    {GR_METHOD_SET_STR,         (gr_funcptr) _gr_fmpz_poly_set_str}, */
    {GR_METHOD_GET_UI,          (gr_funcptr) _gr_fmpz_poly_get_ui},
    {GR_METHOD_GET_SI,          (gr_funcptr) _gr_fmpz_poly_get_si},
    {GR_METHOD_GET_FMPZ,        (gr_funcptr) _gr_fmpz_poly_get_fmpz},
    {GR_METHOD_GET_FMPQ,        (gr_funcptr) _gr_fmpz_poly_get_fmpq},
    {GR_METHOD_NEG,             (gr_funcptr) _gr_fmpz_poly_neg},
    {GR_METHOD_ADD,             (gr_funcptr) _gr_fmpz_poly_add},
    {GR_METHOD_SUB,             (gr_funcptr) _gr_fmpz_poly_sub},
    {GR_METHOD_MUL,             (gr_funcptr) _gr_fmpz_poly_mul},
    {GR_METHOD_MUL_UI,          (gr_funcptr) _gr_fmpz_poly_mul_ui},
    {GR_METHOD_MUL_SI,          (gr_funcptr) _gr_fmpz_poly_mul_si},
    {GR_METHOD_MUL_FMPZ,        (gr_funcptr) _gr_fmpz_poly_mul_fmpz},
    {GR_METHOD_ADDMUL_UI,       (gr_funcptr) _gr_fmpz_poly_addmul_ui},
    {GR_METHOD_ADDMUL_SI,       (gr_funcptr) _gr_fmpz_poly_addmul_si},
    {GR_METHOD_ADDMUL_FMPZ,     (gr_funcptr) _gr_fmpz_poly_addmul_fmpz},
/*    {GR_METHOD_SUBMUL_UI,       (gr_funcptr) _gr_fmpz_poly_submul_ui}, */
/*    {GR_METHOD_SUBMUL_SI,       (gr_funcptr) _gr_fmpz_poly_submul_si}, */
    {GR_METHOD_SUBMUL_FMPZ,     (gr_funcptr) _gr_fmpz_poly_submul_fmpz},
    {GR_METHOD_MUL_TWO,         (gr_funcptr) _gr_fmpz_poly_mul_two},
    {GR_METHOD_SQR,             (gr_funcptr) _gr_fmpz_poly_sqr},
    {GR_METHOD_DIV,             (gr_funcptr) _gr_fmpz_poly_div},
    {GR_METHOD_DIVEXACT,        (gr_funcptr) _gr_fmpz_poly_divexact},
    {GR_METHOD_DIVEXACT_UI,     (gr_funcptr) _gr_fmpz_poly_divexact_ui},
    {GR_METHOD_DIVEXACT_SI,     (gr_funcptr) _gr_fmpz_poly_divexact_si},
    {GR_METHOD_DIVEXACT_FMPZ,   (gr_funcptr) _gr_fmpz_poly_divexact_fmpz},
    {GR_METHOD_IS_INVERTIBLE,   (gr_funcptr) _gr_fmpz_poly_is_invertible},
    {GR_METHOD_INV,             (gr_funcptr) _gr_fmpz_poly_inv},
    {GR_METHOD_DIVIDES,         (gr_funcptr) _gr_fmpz_poly_divides},
    {GR_METHOD_EUCLIDEAN_DIV,   (gr_funcptr) _gr_fmpz_poly_euclidean_div},
    {GR_METHOD_EUCLIDEAN_REM,   (gr_funcptr) _gr_fmpz_poly_euclidean_rem},
    {GR_METHOD_EUCLIDEAN_DIVREM,(gr_funcptr) _gr_fmpz_poly_euclidean_divrem},
    {GR_METHOD_GCD,             (gr_funcptr) _gr_fmpz_poly_gcd},
    {GR_METHOD_LCM,             (gr_funcptr) _gr_fmpz_poly_lcm},
    {GR_METHOD_POW_UI,          (gr_funcptr) _gr_fmpz_poly_pow_ui},
    {GR_METHOD_POW_SI,          (gr_funcptr) _gr_fmpz_poly_pow_si},
    {GR_METHOD_POW_FMPZ,        (gr_funcptr) _gr_fmpz_poly_pow_fmpz},
    {GR_METHOD_IS_SQUARE,       (gr_funcptr) _gr_fmpz_poly_is_square},
    {GR_METHOD_SQRT,            (gr_funcptr) _gr_fmpz_poly_sqrt},
    {GR_METHOD_RSQRT,           (gr_funcptr) _gr_fmpz_poly_rsqrt},
    {GR_METHOD_FACTOR,          (gr_funcptr) _gr_fmpz_poly_factor},
    {0,                         (gr_funcptr) NULL},
};

void
gr_ctx_init_fmpz_poly(gr_ctx_t ctx)
{
    ctx->which_ring = GR_CTX_FMPZ_POLY;
    ctx->sizeof_elem = sizeof(fmpz_poly_struct);
    ctx->size_limit = WORD_MAX;
    ctx->methods = _fmpz_poly_methods;

    POLYNOMIAL_CTX(ctx)->base_ring = NULL;
    POLYNOMIAL_CTX(ctx)->degree_limit = WORD_MAX;
    POLYNOMIAL_CTX(ctx)->var = (char *) default_var;

    if (!_fmpz_poly_methods_initialized)
    {
        gr_method_tab_init(_fmpz_poly_methods, _fmpz_poly_methods_input);
        _fmpz_poly_methods_initialized = 1;
    }
}
