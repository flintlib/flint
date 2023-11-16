/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Polynomials over generic rings */

#include <stdlib.h>
#include <string.h>
#include "fmpz.h"
#include "ulong_extras.h"
#include "gr.h"
#include "gr_poly.h"

#ifdef __GNUC__
# define strcmp __builtin_strcmp
#else
# include <string.h>
#endif

static const char * default_var = "x";

void
polynomial_init(gr_poly_t res, gr_ctx_t ctx)
{
    gr_poly_init(res, POLYNOMIAL_ELEM_CTX(ctx));
}

int polynomial_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    gr_stream_write(out, "Ring of polynomials over ");
    gr_ctx_write(out, POLYNOMIAL_ELEM_CTX(ctx));
    return GR_SUCCESS;
}

int _gr_gr_poly_ctx_set_gen_name(gr_ctx_t ctx, const char * s)
{
    slong len;
    len = strlen(s);

    if (POLYNOMIAL_CTX(ctx)->var == default_var)
        POLYNOMIAL_CTX(ctx)->var = NULL;

    POLYNOMIAL_CTX(ctx)->var = flint_realloc(POLYNOMIAL_CTX(ctx)->var, len + 1);
    memcpy(POLYNOMIAL_CTX(ctx)->var, s, len + 1);
    return GR_SUCCESS;
}

void
polynomial_ctx_clear(gr_ctx_t ctx)
{
    if (POLYNOMIAL_CTX(ctx)->var != default_var)
    {
        flint_free(POLYNOMIAL_CTX(ctx)->var);
    }
}

truth_t
polynomial_ctx_is_ring(gr_ctx_t ctx)
{
    return gr_ctx_is_ring(POLYNOMIAL_ELEM_CTX(ctx));
}

truth_t
polynomial_ctx_is_commutative_ring(gr_ctx_t ctx)
{
    return gr_ctx_is_commutative_ring(POLYNOMIAL_ELEM_CTX(ctx));
}

truth_t
polynomial_ctx_is_integral_domain(gr_ctx_t ctx)
{
    return gr_ctx_is_integral_domain(POLYNOMIAL_ELEM_CTX(ctx));
}

truth_t
polynomial_ctx_is_threadsafe(gr_ctx_t ctx)
{
    return gr_ctx_is_threadsafe(POLYNOMIAL_ELEM_CTX(ctx));
}


void
polynomial_clear(gr_poly_t res, gr_ctx_t ctx)
{
    gr_poly_clear(res, POLYNOMIAL_ELEM_CTX(ctx));
}

void
polynomial_swap(gr_poly_t poly1, gr_poly_t poly2, gr_ctx_t ctx)
{
    gr_poly_swap(poly1, poly2, POLYNOMIAL_ELEM_CTX(ctx));
}

void
polynomial_set_shallow(gr_poly_t res, const gr_poly_t x, const gr_ctx_t ctx)
{
    *res = *x;
}

int
polynomial_write(gr_stream_t out, gr_poly_t poly, gr_ctx_t ctx)
{
    /* todo */
    if (poly->length == 0)
    {
        gr_stream_write(out, "0");
        return GR_SUCCESS;
    }

    return gr_poly_write(out, poly, POLYNOMIAL_CTX(ctx)->var, POLYNOMIAL_ELEM_CTX(ctx));
}

int
polynomial_randtest(gr_poly_t res, flint_rand_t state, gr_ctx_t ctx)
{
    return gr_poly_randtest(res, state, n_randint(state, 5), POLYNOMIAL_ELEM_CTX(ctx));
}

truth_t
polynomial_equal(const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx)
{
    return gr_poly_equal(poly1, poly2, POLYNOMIAL_ELEM_CTX(ctx));
}

int
polynomial_set(gr_poly_t res, const gr_poly_t mat, gr_ctx_t ctx)
{
    return gr_poly_set(res, mat, POLYNOMIAL_ELEM_CTX(ctx));
}

int
polynomial_set_si(gr_poly_t res, slong v, gr_ctx_t ctx)
{
    return gr_poly_set_si(res, v, POLYNOMIAL_ELEM_CTX(ctx));
}

int
polynomial_set_ui(gr_poly_t res, ulong v, gr_ctx_t ctx)
{
    return gr_poly_set_ui(res, v, POLYNOMIAL_ELEM_CTX(ctx));
}

int
polynomial_set_fmpz(gr_poly_t res, const fmpz_t v, gr_ctx_t ctx)
{
    return gr_poly_set_fmpz(res, v, POLYNOMIAL_ELEM_CTX(ctx));
}

int
polynomial_set_fmpq(gr_poly_t res, const fmpq_t v, gr_ctx_t ctx)
{
    return gr_poly_set_fmpq(res, v, POLYNOMIAL_ELEM_CTX(ctx));
}

#include "fmpz_poly.h"

int
polynomial_set_other(gr_poly_t res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
{
    if (x_ctx == ctx)
    {
        return polynomial_set(res, x, ctx);
    }
    else if (x_ctx == POLYNOMIAL_ELEM_CTX(ctx))
    {
        return gr_poly_set_scalar(res, x, x_ctx);
    }
    else if (x_ctx->which_ring == GR_CTX_GR_POLY && !strcmp(POLYNOMIAL_CTX(x_ctx)->var, POLYNOMIAL_CTX(ctx)->var))
    {
        return gr_poly_set_gr_poly_other(res, x, POLYNOMIAL_ELEM_CTX(x_ctx), POLYNOMIAL_ELEM_CTX(ctx));
    }
    else if (x_ctx->which_ring == GR_CTX_FMPZ_POLY)
    {
        return gr_poly_set_fmpz_poly(res, x, POLYNOMIAL_ELEM_CTX(ctx));
    }
    else if (x_ctx->which_ring == GR_CTX_FMPQ_POLY)
    {
        return gr_poly_set_fmpq_poly(res, x, POLYNOMIAL_ELEM_CTX(ctx));
    }
    else if (x_ctx->which_ring == GR_CTX_GR_VEC)
    {
        gr_poly_t tmp;
        tmp->coeffs = ((gr_vec_struct *) x)->entries;
        tmp->length = ((gr_vec_struct *) x)->length;

        return gr_poly_set_gr_poly_other(res, tmp, VECTOR_CTX(x_ctx)->base_ring, POLYNOMIAL_ELEM_CTX(ctx));
    }
    else
    {
        int status = GR_SUCCESS;

        gr_poly_fit_length(res, 1, POLYNOMIAL_ELEM_CTX(ctx));
        status = gr_set_other(res->coeffs, x, x_ctx, POLYNOMIAL_ELEM_CTX(ctx));
        if (status == GR_SUCCESS)
        {
            _gr_poly_set_length(res, 1, POLYNOMIAL_ELEM_CTX(ctx));
            _gr_poly_normalise(res, POLYNOMIAL_ELEM_CTX(ctx));
        }
        else
            _gr_poly_set_length(res, 0, POLYNOMIAL_ELEM_CTX(ctx));
        return status;
    }
}

int
polynomial_zero(gr_poly_t res, gr_ctx_t ctx)
{
    return gr_poly_zero(res, POLYNOMIAL_ELEM_CTX(ctx));
}

int
polynomial_one(gr_poly_t res, gr_ctx_t ctx)
{
    return gr_poly_one(res, POLYNOMIAL_ELEM_CTX(ctx));
}

int
polynomial_neg_one(gr_poly_t res, gr_ctx_t ctx)
{
    return gr_poly_neg_one(res, POLYNOMIAL_ELEM_CTX(ctx));
}

int
polynomial_gen(gr_poly_t res, gr_ctx_t ctx)
{
    return gr_poly_gen(res, POLYNOMIAL_ELEM_CTX(ctx));
}

/*
truth_t
polynomial_is_zero(const gr_poly_t poly, gr_ctx_t ctx)
{
    return gr_poly_is_zero(poly, POLYNOMIAL_ELEM_CTX(ctx));
}

truth_t
polynomial_is_one(const gr_poly_t poly, gr_ctx_t ctx)
{
    return gr_poly_is_one(poly, POLYNOMIAL_ELEM_CTX(ctx));
}

truth_t
polynomial_is_neg_one(const gr_poly_t poly, gr_ctx_t ctx)
{
    return gr_poly_is_neg_one(poly, POLYNOMIAL_ELEM_CTX(ctx));
}
*/

int
polynomial_neg(gr_poly_t res, const gr_poly_t mat, gr_ctx_t ctx)
{
    return gr_poly_neg(res, mat, POLYNOMIAL_ELEM_CTX(ctx));
}

int
polynomial_add(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx)
{
    return gr_poly_add(res, poly1, poly2, POLYNOMIAL_ELEM_CTX(ctx));
}

int
polynomial_sub(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx)
{
    return gr_poly_sub(res, poly1, poly2, POLYNOMIAL_ELEM_CTX(ctx));
}

int
polynomial_mul(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx)
{
    if (POLYNOMIAL_CTX(ctx)->degree_limit != WORD_MAX)
    {
        if (poly1->length != 0 && poly2->length != 0 &&
            poly1->length + poly2->length > POLYNOMIAL_CTX(ctx)->degree_limit)
            return GR_UNABLE;
    }

    return gr_poly_mul(res, poly1, poly2, POLYNOMIAL_ELEM_CTX(ctx));
}

/* todo */
int
polynomial_div(gr_poly_t res, const gr_poly_t x, const gr_poly_t y, const gr_ctx_t ctx)
{
    gr_poly_t r;
    int status;
    gr_poly_init(r, POLYNOMIAL_ELEM_CTX(ctx));
    status = gr_poly_divrem(res, r, x, y, POLYNOMIAL_ELEM_CTX(ctx));

    if (status == GR_SUCCESS)
    {
        truth_t is_zero = gr_poly_is_zero(r, POLYNOMIAL_ELEM_CTX(ctx));

        if (is_zero == T_FALSE)
            status = GR_DOMAIN;
        if (is_zero == T_UNKNOWN)
            status = GR_UNABLE;
    }

    gr_poly_clear(r, POLYNOMIAL_ELEM_CTX(ctx));
    return status;
}

int
polynomial_euclidean_div(gr_poly_t res, const gr_poly_t x, const gr_poly_t y, const gr_ctx_t ctx)
{
    gr_poly_t r;
    int status;
    gr_poly_init(r, POLYNOMIAL_ELEM_CTX(ctx));
    status = gr_poly_divrem(res, r, x, y, POLYNOMIAL_ELEM_CTX(ctx));
    gr_poly_clear(r, POLYNOMIAL_ELEM_CTX(ctx));
    return status;
}

int
polynomial_euclidean_rem(gr_poly_t res, const gr_poly_t x, const gr_poly_t y, const gr_ctx_t ctx)
{
    gr_poly_t q;
    int status;
    gr_poly_init(q, POLYNOMIAL_ELEM_CTX(ctx));
    status = gr_poly_divrem(q, res, x, y, POLYNOMIAL_ELEM_CTX(ctx));
    gr_poly_clear(q, POLYNOMIAL_ELEM_CTX(ctx));
    return status;
}

int
polynomial_euclidean_divrem(gr_poly_t res1, gr_poly_t res2, const gr_poly_t x, const gr_poly_t y, const gr_ctx_t ctx)
{
    return gr_poly_divrem(res1, res2, x, y, POLYNOMIAL_ELEM_CTX(ctx));
}

int
polynomial_inv(gr_poly_t res, const gr_poly_t poly, gr_ctx_t ctx)
{
    return gr_poly_inv(res, poly, POLYNOMIAL_ELEM_CTX(ctx));
}

int
polynomial_pow_ui(gr_poly_t res, const gr_poly_t poly, ulong exp, gr_ctx_t ctx)
{
    return gr_poly_pow_ui(res, poly, exp, POLYNOMIAL_ELEM_CTX(ctx));
}

int
polynomial_pow_fmpz(gr_poly_t res, const gr_poly_t poly, const fmpz_t exp, gr_ctx_t ctx)
{
    return gr_poly_pow_fmpz(res, poly, exp, POLYNOMIAL_ELEM_CTX(ctx));
}

int
polynomial_pow_si(gr_poly_t res, const gr_poly_t poly, slong exp, gr_ctx_t ctx)
{
    int status;
    fmpz_t t;
    fmpz_init_set_si(t, exp);
    status = gr_poly_pow_fmpz(res, poly, t, POLYNOMIAL_ELEM_CTX(ctx));
    fmpz_clear(t);
    return status;
}

int
polynomial_gcd(gr_poly_t res, const gr_poly_t x, const gr_poly_t y, const gr_ctx_t ctx)
{
    return gr_poly_gcd(res, x, y, POLYNOMIAL_ELEM_CTX(ctx));
}


int _gr_poly_methods_initialized = 0;

gr_static_method_table _gr_poly_methods;

gr_method_tab_input _gr_poly_methods_input[] =
{
    {GR_METHOD_CTX_WRITE,   (gr_funcptr) polynomial_ctx_write},
    {GR_METHOD_CTX_CLEAR,   (gr_funcptr) polynomial_ctx_clear},

    {GR_METHOD_CTX_IS_RING,     (gr_funcptr) polynomial_ctx_is_ring},
    {GR_METHOD_CTX_IS_COMMUTATIVE_RING, (gr_funcptr) polynomial_ctx_is_commutative_ring},
    {GR_METHOD_CTX_IS_INTEGRAL_DOMAIN,  (gr_funcptr) polynomial_ctx_is_integral_domain},
    {GR_METHOD_CTX_IS_FIELD,            (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_THREADSAFE,       (gr_funcptr) polynomial_ctx_is_threadsafe},
    {GR_METHOD_CTX_SET_GEN_NAME, (gr_funcptr) _gr_gr_poly_ctx_set_gen_name},

    {GR_METHOD_INIT,        (gr_funcptr) polynomial_init},
    {GR_METHOD_CLEAR,       (gr_funcptr) polynomial_clear},
    {GR_METHOD_SWAP,        (gr_funcptr) polynomial_swap},
    {GR_METHOD_SET_SHALLOW, (gr_funcptr) polynomial_set_shallow},
    {GR_METHOD_RANDTEST,    (gr_funcptr) polynomial_randtest},
    {GR_METHOD_WRITE,       (gr_funcptr) polynomial_write},
    {GR_METHOD_ZERO,        (gr_funcptr) polynomial_zero},
    {GR_METHOD_ONE,         (gr_funcptr) polynomial_one},
    {GR_METHOD_NEG_ONE,     (gr_funcptr) polynomial_neg_one},

    {GR_METHOD_GEN,         (gr_funcptr) polynomial_gen},

/*
    {GR_METHOD_IS_ZERO,     (gr_funcptr) polynomial_is_zero},
    {GR_METHOD_IS_ONE,      (gr_funcptr) polynomial_is_one},
    {GR_METHOD_IS_NEG_ONE,  (gr_funcptr) polynomial_is_neg_one},
*/
    {GR_METHOD_EQUAL,       (gr_funcptr) polynomial_equal},
    {GR_METHOD_SET,         (gr_funcptr) polynomial_set},
    {GR_METHOD_SET_UI,      (gr_funcptr) polynomial_set_ui},
    {GR_METHOD_SET_SI,      (gr_funcptr) polynomial_set_si},
    {GR_METHOD_SET_FMPZ,    (gr_funcptr) polynomial_set_fmpz},
    {GR_METHOD_SET_FMPQ,    (gr_funcptr) polynomial_set_fmpq},
    {GR_METHOD_SET_OTHER,   (gr_funcptr) polynomial_set_other},
    {GR_METHOD_NEG,         (gr_funcptr) polynomial_neg},
    {GR_METHOD_ADD,         (gr_funcptr) polynomial_add},
    {GR_METHOD_SUB,         (gr_funcptr) polynomial_sub},
    {GR_METHOD_MUL,         (gr_funcptr) polynomial_mul},
    {GR_METHOD_POW_UI,      (gr_funcptr) polynomial_pow_ui},
    {GR_METHOD_POW_SI,      (gr_funcptr) polynomial_pow_si},
    {GR_METHOD_POW_FMPZ,    (gr_funcptr) polynomial_pow_fmpz},
    {GR_METHOD_DIV,         (gr_funcptr) polynomial_div},
    {GR_METHOD_INV,         (gr_funcptr) polynomial_inv},

    {GR_METHOD_EUCLIDEAN_DIV,         (gr_funcptr) polynomial_euclidean_div},
    {GR_METHOD_EUCLIDEAN_REM,         (gr_funcptr) polynomial_euclidean_rem},
    {GR_METHOD_EUCLIDEAN_DIVREM,      (gr_funcptr) polynomial_euclidean_divrem},

    {GR_METHOD_GCD,         (gr_funcptr) polynomial_gcd},

    {0,                     (gr_funcptr) NULL},
};

void
gr_ctx_init_gr_poly(gr_ctx_t ctx, gr_ctx_t base_ring)
{
    ctx->which_ring = GR_CTX_GR_POLY;
    ctx->sizeof_elem = sizeof(gr_poly_struct);
    ctx->size_limit = WORD_MAX;

    POLYNOMIAL_CTX(ctx)->base_ring = (gr_ctx_struct *) base_ring;
    POLYNOMIAL_CTX(ctx)->degree_limit = WORD_MAX;
    POLYNOMIAL_CTX(ctx)->var = (char *) default_var;

    ctx->methods = _gr_poly_methods;

    if (!_gr_poly_methods_initialized)
    {
        gr_method_tab_init(_gr_poly_methods, _gr_poly_methods_input);
        _gr_poly_methods_initialized = 1;
    }
}
