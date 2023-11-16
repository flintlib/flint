/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "gr.h"
#include "gr_vec.h"
#include "gr_generic.h"
#include "fmpz_mpoly.h"
#include "fmpz_mpoly.h"
#include "fmpz_mpoly_factor.h"

typedef struct
{
    fmpz_mpoly_ctx_t mctx;
}
_gr_fmpz_mpoly_ctx_t;

#define MPOLYNOMIAL_CTX(ring_ctx) ((_gr_fmpz_mpoly_ctx_t *)(GR_CTX_DATA_AS_PTR(ring_ctx)))
#define MPOLYNOMIAL_MCTX(ring_ctx) (MPOLYNOMIAL_CTX(ring_ctx)->mctx)

int _gr_fmpz_mpoly_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    gr_stream_write(out, "Ring of multivariate polynomials over Integer ring (fmpz)");
    gr_stream_write(out, " in ");
    gr_stream_write_si(out, MPOLYNOMIAL_MCTX(ctx)->minfo->nvars);
    gr_stream_write(out, " variables");
    if (MPOLYNOMIAL_MCTX(ctx)->minfo->ord == ORD_LEX)
        gr_stream_write(out, ", lex order");
    else if (MPOLYNOMIAL_MCTX(ctx)->minfo->ord == ORD_DEGLEX)
        gr_stream_write(out, ", deglex order");
    else if (MPOLYNOMIAL_MCTX(ctx)->minfo->ord == ORD_DEGREVLEX)
        gr_stream_write(out, ", degrevlex order");
    return GR_SUCCESS;
}

void
_gr_fmpz_mpoly_ctx_clear(gr_ctx_t ctx)
{
    fmpz_mpoly_ctx_clear(MPOLYNOMIAL_MCTX(ctx));
    flint_free(GR_CTX_DATA_AS_PTR(ctx));
}

void
_gr_fmpz_mpoly_init(fmpz_mpoly_t res, gr_ctx_t ctx)
{
    fmpz_mpoly_init(res, MPOLYNOMIAL_MCTX(ctx));
}

void
_gr_fmpz_mpoly_clear(fmpz_mpoly_t res, gr_ctx_t ctx)
{
    fmpz_mpoly_clear(res, MPOLYNOMIAL_MCTX(ctx));
}

void
_gr_fmpz_mpoly_swap(fmpz_mpoly_t poly1, fmpz_mpoly_t poly2, gr_ctx_t ctx)
{
    fmpz_mpoly_swap(poly1, poly2, MPOLYNOMIAL_MCTX(ctx));
}

void
_gr_fmpz_mpoly_set_shallow(fmpz_mpoly_t res, const fmpz_mpoly_t poly, gr_ctx_t ctx)
{
    *res = *poly;
}

int
_gr_fmpz_mpoly_randtest(fmpz_mpoly_t res, flint_rand_t state, gr_ctx_t ctx)
{
    slong bits;

    if (n_randint(state, 10) != 0)
        bits = 10;
    else
        bits = 100;

    fmpz_mpoly_randtest_bits(res, state, n_randint(state, 5), bits, 1 + n_randint(state, 3), MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mpoly_randtest_small(fmpz_mpoly_t res, flint_rand_t state, gr_ctx_t ctx)
{
    fmpz_mpoly_randtest_bits(res, state, n_randint(state, 3), 3, 1 + n_randint(state, 3), MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mpoly_write(gr_stream_t out, fmpz_mpoly_t poly, gr_ctx_t ctx)
{
/*
    if (out->fp != NULL)
        fmpz_mpoly_fprint_pretty(out->fp, poly, NULL, MPOLYNOMIAL_MCTX(ctx));
    else
*/
    gr_stream_write_free(out, fmpz_mpoly_get_str_pretty(poly, NULL, MPOLYNOMIAL_MCTX(ctx)));

    /* todo: error handling */
    return GR_SUCCESS;
}

truth_t
_gr_fmpz_mpoly_equal(const fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, gr_ctx_t ctx)
{
    return fmpz_mpoly_equal(poly1, poly2, MPOLYNOMIAL_MCTX(ctx)) ? T_TRUE : T_FALSE;
}

truth_t
_gr_fmpz_mpoly_is_zero(const fmpz_mpoly_t poly, gr_ctx_t ctx)
{
    return fmpz_mpoly_is_zero(poly, MPOLYNOMIAL_MCTX(ctx)) ? T_TRUE : T_FALSE;
}

truth_t
_gr_fmpz_mpoly_is_one(const fmpz_mpoly_t poly, gr_ctx_t ctx)
{
    return fmpz_mpoly_is_one(poly, MPOLYNOMIAL_MCTX(ctx)) ? T_TRUE : T_FALSE;
}

int
_gr_fmpz_mpoly_zero(fmpz_mpoly_t res, gr_ctx_t ctx)
{
    fmpz_mpoly_zero(res, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mpoly_one(fmpz_mpoly_t res, gr_ctx_t ctx)
{
    fmpz_mpoly_one(res, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

truth_t
_gr_fmpz_mpoly_gens(gr_vec_t res, gr_ctx_t ctx)
{
    slong i, n;

    n = MPOLYNOMIAL_MCTX(ctx)->minfo->nvars;

    gr_vec_set_length(res, n, ctx);
    for (i = 0; i < n; i++)
        fmpz_mpoly_gen(((fmpz_mpoly_struct *) res->entries) + i, i, MPOLYNOMIAL_MCTX(ctx));

    return GR_SUCCESS;
}

int
_gr_fmpz_mpoly_set(fmpz_mpoly_t res, const fmpz_mpoly_t mat, gr_ctx_t ctx)
{
    fmpz_mpoly_set(res, mat, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mpoly_set_si(fmpz_mpoly_t res, slong v, gr_ctx_t ctx)
{
    fmpz_mpoly_set_si(res, v, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mpoly_set_ui(fmpz_mpoly_t res, ulong v, gr_ctx_t ctx)
{
    fmpz_mpoly_set_ui(res, v, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mpoly_set_fmpz(fmpz_mpoly_t res, const fmpz_t v, gr_ctx_t ctx)
{
    fmpz_mpoly_set_fmpz(res, v, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mpoly_set_fmpq(fmpz_mpoly_t res, const fmpq_t v, gr_ctx_t ctx)
{
    if (!fmpz_is_one(fmpq_denref(v)))
        return GR_DOMAIN;

    fmpz_mpoly_set_fmpz(res, fmpq_numref(v), MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mpoly_set_other(fmpz_mpoly_t res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
{
    if (x_ctx->which_ring == GR_CTX_FMPZ_MPOLY)
    {
        /* fmpz_mpoly_set_fmpz_poly */
        if (MPOLYNOMIAL_MCTX(ctx)->minfo->nvars == MPOLYNOMIAL_MCTX(x_ctx)->minfo->nvars &&
            MPOLYNOMIAL_MCTX(ctx)->minfo->ord == MPOLYNOMIAL_MCTX(x_ctx)->minfo->ord)
        {
            fmpz_mpoly_set(res, x, MPOLYNOMIAL_MCTX(ctx));
            return GR_SUCCESS;
        }
    }

    return gr_generic_set_other(res, x, x_ctx, ctx);
}

int
_gr_fmpz_mpoly_neg(fmpz_mpoly_t res, const fmpz_mpoly_t mat, gr_ctx_t ctx)
{
    fmpz_mpoly_neg(res, mat, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mpoly_add(fmpz_mpoly_t res, const fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, gr_ctx_t ctx)
{
    if (poly1->length + poly2->length > ctx->size_limit)
    {
        fmpz_mpoly_zero(res, MPOLYNOMIAL_MCTX(ctx));
        return GR_UNABLE;
    }

    fmpz_mpoly_add(res, poly1, poly2, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mpoly_add_si(fmpz_mpoly_t res, const fmpz_mpoly_t poly1, slong c, gr_ctx_t ctx)
{
    fmpz_mpoly_add_si(res, poly1, c, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mpoly_add_ui(fmpz_mpoly_t res, const fmpz_mpoly_t poly1, ulong c, gr_ctx_t ctx)
{
    fmpz_mpoly_add_ui(res, poly1, c, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mpoly_add_fmpz(fmpz_mpoly_t res, const fmpz_mpoly_t poly1, const fmpz_t c, gr_ctx_t ctx)
{
    fmpz_mpoly_add_fmpz(res, poly1, c, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mpoly_sub(fmpz_mpoly_t res, const fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, gr_ctx_t ctx)
{
    if (poly1->length + poly2->length > ctx->size_limit)
    {
        fmpz_mpoly_zero(res, MPOLYNOMIAL_MCTX(ctx));
        return GR_UNABLE;
    }

    fmpz_mpoly_sub(res, poly1, poly2, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mpoly_sub_si(fmpz_mpoly_t res, const fmpz_mpoly_t poly1, slong c, gr_ctx_t ctx)
{
    fmpz_mpoly_sub_si(res, poly1, c, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mpoly_sub_ui(fmpz_mpoly_t res, const fmpz_mpoly_t poly1, ulong c, gr_ctx_t ctx)
{
    fmpz_mpoly_sub_ui(res, poly1, c, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mpoly_sub_fmpz(fmpz_mpoly_t res, const fmpz_mpoly_t poly1, const fmpz_t c, gr_ctx_t ctx)
{
    fmpz_mpoly_sub_fmpz(res, poly1, c, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mpoly_mul(fmpz_mpoly_t res, const fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, gr_ctx_t ctx)
{
    if (poly1->length * poly2->length > ctx->size_limit)  /* todo: * can overflow */
    {
        fmpz_mpoly_zero(res, MPOLYNOMIAL_MCTX(ctx));
        return GR_UNABLE;
    }

    fmpz_mpoly_mul(res, poly1, poly2, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mpoly_mul_si(fmpz_mpoly_t res, const fmpz_mpoly_t poly1, slong c, gr_ctx_t ctx)
{
    fmpz_mpoly_scalar_mul_si(res, poly1, c, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mpoly_mul_ui(fmpz_mpoly_t res, const fmpz_mpoly_t poly1, ulong c, gr_ctx_t ctx)
{
    fmpz_mpoly_scalar_mul_ui(res, poly1, c, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mpoly_mul_fmpz(fmpz_mpoly_t res, const fmpz_mpoly_t poly1, const fmpz_t c, gr_ctx_t ctx)
{
    fmpz_mpoly_scalar_mul_fmpz(res, poly1, c, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mpoly_div(fmpz_mpoly_t res, const fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, gr_ctx_t ctx)
{
    if (poly2->length == 0)
        return GR_DOMAIN;

    if (fmpz_mpoly_divides(res, poly1, poly2, MPOLYNOMIAL_MCTX(ctx)))
        return GR_SUCCESS;
    else
        return GR_DOMAIN;
}

int
_gr_fmpz_mpoly_divexact(fmpz_mpoly_t res, const fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, gr_ctx_t ctx)
{
    if (poly2->length == 0)
        return GR_DOMAIN;

    fmpz_mpoly_div(res, poly1, poly2, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mpoly_divexact_si(fmpz_mpoly_t res, const fmpz_mpoly_t poly1, slong c, gr_ctx_t ctx)
{
    if (c == 0)
        return GR_DOMAIN;

    fmpz_mpoly_scalar_divexact_si(res, poly1, c, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mpoly_divexact_ui(fmpz_mpoly_t res, const fmpz_mpoly_t poly1, ulong c, gr_ctx_t ctx)
{
    if (c == 0)
        return GR_DOMAIN;

    fmpz_mpoly_scalar_divexact_ui(res, poly1, c, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mpoly_divexact_fmpz(fmpz_mpoly_t res, const fmpz_mpoly_t poly1, const fmpz_t c, gr_ctx_t ctx)
{
    if (fmpz_is_zero(c))
        return GR_DOMAIN;

    fmpz_mpoly_scalar_divexact_fmpz(res, poly1, c, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

truth_t
_gr_fmpz_mpoly_is_invertible(const fmpz_mpoly_t c, gr_ctx_t ctx)
{
    if (c->length == 1 && fmpz_mpoly_is_fmpz(c, MPOLYNOMIAL_MCTX(ctx)) && fmpz_is_pm1(c->coeffs))
        return T_TRUE;

    return T_FALSE;
}

int
_gr_fmpz_mpoly_inv(fmpz_mpoly_t res, const fmpz_mpoly_t c, gr_ctx_t ctx)
{
    if (c->length == 1 && fmpz_mpoly_is_fmpz(c, MPOLYNOMIAL_MCTX(ctx)) && fmpz_is_pm1(c->coeffs))
    {
        fmpz_mpoly_set(res, c, MPOLYNOMIAL_MCTX(ctx));
        return GR_SUCCESS;
    }

    return GR_DOMAIN;
}

int
_gr_fmpz_mpoly_pow_ui(fmpz_mpoly_t res, const fmpz_mpoly_t poly1, ulong c, gr_ctx_t ctx)
{
    /* todo: size limit */

    if (fmpz_mpoly_pow_ui(res, poly1, c, MPOLYNOMIAL_MCTX(ctx)))
        return GR_SUCCESS;
    else
        return GR_UNABLE;
}

int
_gr_fmpz_mpoly_pow_fmpz(fmpz_mpoly_t res, const fmpz_mpoly_t poly1, const fmpz_t c, gr_ctx_t ctx)
{
    /* todo: size limit */

    if (fmpz_sgn(c) < 0)
    {
        int status;

        status = gr_inv(res, poly1, ctx);

        if (status == GR_SUCCESS)
        {
            fmpz_t e;
            fmpz_init(e);
            fmpz_neg(e, c);
            status = _gr_fmpz_mpoly_pow_fmpz(res, res, e, ctx);
            fmpz_clear(e);
        }

        return status;
    }

    if (fmpz_mpoly_pow_fmpz(res, poly1, c, MPOLYNOMIAL_MCTX(ctx)))
        return GR_SUCCESS;
    else
        return GR_UNABLE;
}

truth_t
_gr_fmpz_mpoly_is_square(const fmpz_mpoly_t poly, gr_ctx_t ctx)
{
    return fmpz_mpoly_is_square(poly, MPOLYNOMIAL_MCTX(ctx)) ? T_TRUE : T_FALSE;
}

int
_gr_fmpz_mpoly_sqrt(fmpz_mpoly_t res, const fmpz_mpoly_t poly, gr_ctx_t ctx)
{
    if (fmpz_mpoly_sqrt(res, poly, MPOLYNOMIAL_MCTX(ctx)))
        return GR_SUCCESS;
    else
        return GR_DOMAIN;
}

int
_gr_fmpz_mpoly_gcd(fmpz_mpoly_t res, const fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, gr_ctx_t ctx)
{
    if (fmpz_mpoly_gcd(res, poly1, poly2, MPOLYNOMIAL_MCTX(ctx)))
        return GR_SUCCESS;

    return GR_DOMAIN;
}

int
_gr_fmpz_mpoly_factor(fmpz_mpoly_t c, gr_vec_t factors, gr_vec_t exponents, gr_srcptr x, int flags, gr_ctx_t ctx)
{
    fmpz_mpoly_factor_t fac;
    gr_ctx_t ZZ;
    slong i;
    int status = GR_SUCCESS;

    fmpz_mpoly_factor_init(fac, MPOLYNOMIAL_MCTX(ctx));

    if (fmpz_mpoly_factor(fac, x, MPOLYNOMIAL_MCTX(ctx)))
    {
        fmpz_mpoly_set_fmpz(c, fac->constant, MPOLYNOMIAL_MCTX(ctx));

        gr_ctx_init_fmpz(ZZ);

        gr_vec_set_length(factors, fac->num, ctx);
        gr_vec_set_length(exponents, fac->num, ZZ);

        for (i = 0; i < fac->num; i++)
        {
            fmpz_mpoly_swap((fmpz_mpoly_struct *) (factors->entries) + i, fac->poly + i, MPOLYNOMIAL_MCTX(ctx));
            fmpz_swap((fmpz *) (exponents->entries) + i, fac->exp + i);
        }

        gr_ctx_clear(ZZ);
    }
    else
    {
        status = GR_UNABLE;
    }

    fmpz_mpoly_factor_clear(fac, MPOLYNOMIAL_MCTX(ctx));

    return status;
}


int _gr_fmpz_mpoly_methods_initialized = 0;

gr_static_method_table _gr_fmpz_mpoly_methods;

gr_method_tab_input _gr_fmpz_mpoly_methods_input[] =
{
    {GR_METHOD_CTX_WRITE,   (gr_funcptr) _gr_fmpz_mpoly_ctx_write},
    {GR_METHOD_CTX_CLEAR,   (gr_funcptr) _gr_fmpz_mpoly_ctx_clear},
    {GR_METHOD_CTX_IS_RING,                     (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_COMMUTATIVE_RING,         (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_INTEGRAL_DOMAIN,          (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_FIELD,                    (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_FINITE,                   (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_FINITE_CHARACTERISTIC,    (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_THREADSAFE,               (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_INIT,        (gr_funcptr) _gr_fmpz_mpoly_init},
    {GR_METHOD_CLEAR,       (gr_funcptr) _gr_fmpz_mpoly_clear},
    {GR_METHOD_SWAP,        (gr_funcptr) _gr_fmpz_mpoly_swap},
    {GR_METHOD_SET_SHALLOW, (gr_funcptr) _gr_fmpz_mpoly_set_shallow},
    {GR_METHOD_RANDTEST,    (gr_funcptr) _gr_fmpz_mpoly_randtest},
    {GR_METHOD_RANDTEST_SMALL,    (gr_funcptr) _gr_fmpz_mpoly_randtest_small},
    {GR_METHOD_WRITE,       (gr_funcptr) _gr_fmpz_mpoly_write},
    {GR_METHOD_ZERO,        (gr_funcptr) _gr_fmpz_mpoly_zero},
    {GR_METHOD_ONE,         (gr_funcptr) _gr_fmpz_mpoly_one},
    {GR_METHOD_IS_ZERO,     (gr_funcptr) _gr_fmpz_mpoly_is_zero},
    {GR_METHOD_IS_ONE,      (gr_funcptr) _gr_fmpz_mpoly_is_one},
    {GR_METHOD_GENS,        (gr_funcptr) _gr_fmpz_mpoly_gens},
    {GR_METHOD_EQUAL,       (gr_funcptr) _gr_fmpz_mpoly_equal},
    {GR_METHOD_SET,         (gr_funcptr) _gr_fmpz_mpoly_set},
    {GR_METHOD_SET_UI,      (gr_funcptr) _gr_fmpz_mpoly_set_ui},
    {GR_METHOD_SET_SI,      (gr_funcptr) _gr_fmpz_mpoly_set_si},
    {GR_METHOD_SET_FMPZ,    (gr_funcptr) _gr_fmpz_mpoly_set_fmpz},
    {GR_METHOD_SET_FMPQ,    (gr_funcptr) _gr_fmpz_mpoly_set_fmpq},
    {GR_METHOD_SET_OTHER,   (gr_funcptr) _gr_fmpz_mpoly_set_other},
    {GR_METHOD_NEG,         (gr_funcptr) _gr_fmpz_mpoly_neg},
    {GR_METHOD_ADD,         (gr_funcptr) _gr_fmpz_mpoly_add},
    {GR_METHOD_ADD_SI,      (gr_funcptr) _gr_fmpz_mpoly_add_si},
    {GR_METHOD_ADD_UI,      (gr_funcptr) _gr_fmpz_mpoly_add_ui},
    {GR_METHOD_ADD_FMPZ,    (gr_funcptr) _gr_fmpz_mpoly_add_fmpz},
    {GR_METHOD_SUB,         (gr_funcptr) _gr_fmpz_mpoly_sub},
    {GR_METHOD_SUB_SI,      (gr_funcptr) _gr_fmpz_mpoly_sub_si},
    {GR_METHOD_SUB_UI,      (gr_funcptr) _gr_fmpz_mpoly_sub_ui},
    {GR_METHOD_SUB_FMPZ,    (gr_funcptr) _gr_fmpz_mpoly_sub_fmpz},
    {GR_METHOD_MUL,         (gr_funcptr) _gr_fmpz_mpoly_mul},
    {GR_METHOD_MUL_SI,      (gr_funcptr) _gr_fmpz_mpoly_mul_si},
    {GR_METHOD_MUL_UI,      (gr_funcptr) _gr_fmpz_mpoly_mul_ui},
    {GR_METHOD_MUL_FMPZ,    (gr_funcptr) _gr_fmpz_mpoly_mul_fmpz},
    {GR_METHOD_DIV,         (gr_funcptr) _gr_fmpz_mpoly_div},
    {GR_METHOD_DIVEXACT,         (gr_funcptr) _gr_fmpz_mpoly_divexact},
    {GR_METHOD_DIVEXACT_SI,      (gr_funcptr) _gr_fmpz_mpoly_divexact_si},
    {GR_METHOD_DIVEXACT_UI,      (gr_funcptr) _gr_fmpz_mpoly_divexact_ui},
    {GR_METHOD_DIVEXACT_FMPZ,    (gr_funcptr) _gr_fmpz_mpoly_divexact_fmpz},
    {GR_METHOD_INV,             (gr_funcptr) _gr_fmpz_mpoly_inv},
    {GR_METHOD_IS_INVERTIBLE,   (gr_funcptr) _gr_fmpz_mpoly_is_invertible},
    {GR_METHOD_POW_UI,      (gr_funcptr) _gr_fmpz_mpoly_pow_ui},
    {GR_METHOD_POW_FMPZ,    (gr_funcptr) _gr_fmpz_mpoly_pow_fmpz},
    {GR_METHOD_SQRT,        (gr_funcptr) _gr_fmpz_mpoly_sqrt},
    {GR_METHOD_IS_SQUARE,   (gr_funcptr) _gr_fmpz_mpoly_is_square},
    {GR_METHOD_GCD,         (gr_funcptr) _gr_fmpz_mpoly_gcd},
    {GR_METHOD_FACTOR,      (gr_funcptr) _gr_fmpz_mpoly_factor},
    {0,                     (gr_funcptr) NULL},
};

void
gr_ctx_init_fmpz_mpoly(gr_ctx_t ctx, slong nvars, const ordering_t ord)
{
    ctx->which_ring = GR_CTX_FMPZ_MPOLY;
    ctx->sizeof_elem = sizeof(fmpz_mpoly_struct);
    GR_CTX_DATA_AS_PTR(ctx) = flint_malloc(sizeof(_gr_fmpz_mpoly_ctx_t));
    ctx->size_limit = WORD_MAX;

    fmpz_mpoly_ctx_init(MPOLYNOMIAL_MCTX(ctx), nvars, ord);

    ctx->methods = _gr_fmpz_mpoly_methods;

    if (!_gr_fmpz_mpoly_methods_initialized)
    {
        gr_method_tab_init(_gr_fmpz_mpoly_methods, _gr_fmpz_mpoly_methods_input);
        _gr_fmpz_mpoly_methods_initialized = 1;
    }
}
