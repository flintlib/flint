/*
    Copyright (C) 2023, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "mpoly.h"
#include "fmpz.h"
#include "fmpq_mpoly_factor.h"
#include "gr.h"
#include "gr/impl.h"
#include "gr_vec.h"
#include "gr_generic.h"

typedef struct
{
    fmpq_mpoly_ctx_t mctx;
    char ** vars;
}
_gr_fmpq_mpoly_ctx_t;

#define MPOLYNOMIAL_CTX(ring_ctx) ((_gr_fmpq_mpoly_ctx_t *)(GR_CTX_DATA_AS_PTR(ring_ctx)))
#define MPOLYNOMIAL_MCTX(ring_ctx) (MPOLYNOMIAL_CTX(ring_ctx)->mctx)

static int _gr_fmpq_mpoly_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    gr_stream_write(out, "Ring of multivariate polynomials over Rational field (fmpq)");
    gr_stream_write(out, " in ");
    gr_stream_write_si(out, MPOLYNOMIAL_MCTX(ctx)->zctx->minfo->nvars);
    gr_stream_write(out, " variables");
    if (MPOLYNOMIAL_MCTX(ctx)->zctx->minfo->ord == ORD_LEX)
        gr_stream_write(out, ", lex order");
    else if (MPOLYNOMIAL_MCTX(ctx)->zctx->minfo->ord == ORD_DEGLEX)
        gr_stream_write(out, ", deglex order");
    else if (MPOLYNOMIAL_MCTX(ctx)->zctx->minfo->ord == ORD_DEGREVLEX)
        gr_stream_write(out, ", degrevlex order");
    return GR_SUCCESS;
}

static void
_gr_fmpq_mpoly_ctx_clear(gr_ctx_t ctx)
{
    if (MPOLYNOMIAL_CTX(ctx)->vars != NULL)
    {
        slong i;
        for (i = 0; i < MPOLYNOMIAL_MCTX(ctx)->zctx->minfo->nvars; i++)
            flint_free(MPOLYNOMIAL_CTX(ctx)->vars[i]);
        flint_free(MPOLYNOMIAL_CTX(ctx)->vars);
    }

    fmpq_mpoly_ctx_clear(MPOLYNOMIAL_MCTX(ctx));
    flint_free(GR_CTX_DATA_AS_PTR(ctx));
}

static int
_gr_fmpq_mpoly_ctx_set_gen_names(gr_ctx_t ctx, const char ** s)
{
    slong i, nvars, len;

    nvars = MPOLYNOMIAL_MCTX(ctx)->zctx->minfo->nvars;

    if (MPOLYNOMIAL_CTX(ctx)->vars == NULL)
    {
        MPOLYNOMIAL_CTX(ctx)->vars = flint_malloc(nvars * sizeof(char *));
        for (i = 0; i < nvars; i++)
            MPOLYNOMIAL_CTX(ctx)->vars[i] = NULL;
    }
    else
    {
        for (i = 0; i < nvars; i++)
            flint_free(MPOLYNOMIAL_CTX(ctx)->vars[i]);
    }

    for (i = 0; i < nvars; i++)
    {
        len = strlen(s[i]);
        MPOLYNOMIAL_CTX(ctx)->vars[i] = flint_realloc(MPOLYNOMIAL_CTX(ctx)->vars[i], len + 1);
        memcpy(MPOLYNOMIAL_CTX(ctx)->vars[i], s[i], len + 1);
    }

    return GR_SUCCESS;
}

static slong
_gr_fmpq_mpoly_ctx_ngens(slong * ngens, gr_ctx_t ctx)
{
    * ngens = MPOLYNOMIAL_MCTX(ctx)->zctx->minfo->nvars;
    return GR_SUCCESS;
}

static int
_gr_fmpq_mpoly_ctx_gen_name(char ** name, slong i, gr_ctx_t ctx)
{
    if (i < 0 || i >= MPOLYNOMIAL_MCTX(ctx)->zctx->minfo->nvars)
        return GR_DOMAIN;

    if (MPOLYNOMIAL_CTX(ctx)->vars == NULL)
        return GR_UNABLE;

    char * var = MPOLYNOMIAL_CTX(ctx)->vars[i];
    size_t len = strlen(var);
    * name = flint_malloc(len + 1);
    if (* name == NULL)
        return GR_UNABLE;
    strncpy(* name, var, len + 1);

    return GR_SUCCESS;
}

static void
_gr_fmpq_mpoly_init(fmpq_mpoly_t res, gr_ctx_t ctx)
{
    fmpq_mpoly_init(res, MPOLYNOMIAL_MCTX(ctx));
}

static void
_gr_fmpq_mpoly_clear(fmpq_mpoly_t res, gr_ctx_t ctx)
{
    fmpq_mpoly_clear(res, MPOLYNOMIAL_MCTX(ctx));
}

static void
_gr_fmpq_mpoly_swap(fmpq_mpoly_t poly1, fmpq_mpoly_t poly2, gr_ctx_t ctx)
{
    fmpq_mpoly_swap(poly1, poly2, MPOLYNOMIAL_MCTX(ctx));
}

static void
_gr_fmpq_mpoly_set_shallow(fmpq_mpoly_t res, const fmpq_mpoly_t poly, gr_ctx_t ctx)
{
    *res = *poly;
}

static int
_gr_fmpq_mpoly_randtest(fmpq_mpoly_t res, flint_rand_t state, gr_ctx_t ctx)
{
    slong bits;

    if (n_randint(state, 10) != 0)
        bits = 10;
    else
        bits = 100;

    fmpq_mpoly_randtest_bits(res, state, n_randint(state, 5), bits, 1 + n_randint(state, 3), MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpq_mpoly_randtest_small(fmpq_mpoly_t res, flint_rand_t state, gr_ctx_t ctx)
{
    fmpq_mpoly_randtest_bits(res, state, n_randint(state, 3), 3, 1 + n_randint(state, 3), MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static slong
_gr_fmpq_mpoly_length(const fmpq_mpoly_t x, gr_ctx_t ctx)
{
    return fmpq_mpoly_length(x, MPOLYNOMIAL_MCTX(ctx));
}

static int
_gr_fmpq_mpoly_write(gr_stream_t out, fmpq_mpoly_t poly, gr_ctx_t ctx)
{
    gr_stream_write_free(out, fmpq_mpoly_get_str_pretty(poly, (const char **) MPOLYNOMIAL_CTX(ctx)->vars, MPOLYNOMIAL_MCTX(ctx)));
    return GR_SUCCESS;
}

static truth_t
_gr_fmpq_mpoly_equal(const fmpq_mpoly_t poly1, const fmpq_mpoly_t poly2, gr_ctx_t ctx)
{
    return fmpq_mpoly_equal(poly1, poly2, MPOLYNOMIAL_MCTX(ctx)) ? T_TRUE : T_FALSE;
}

static truth_t
_gr_fmpq_mpoly_is_zero(const fmpq_mpoly_t poly, gr_ctx_t ctx)
{
    return fmpq_mpoly_is_zero(poly, MPOLYNOMIAL_MCTX(ctx)) ? T_TRUE : T_FALSE;
}

static truth_t
_gr_fmpq_mpoly_is_one(const fmpq_mpoly_t poly, gr_ctx_t ctx)
{
    return fmpq_mpoly_is_one(poly, MPOLYNOMIAL_MCTX(ctx)) ? T_TRUE : T_FALSE;
}

static int
_gr_fmpq_mpoly_zero(fmpq_mpoly_t res, gr_ctx_t ctx)
{
    fmpq_mpoly_zero(res, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpq_mpoly_one(fmpq_mpoly_t res, gr_ctx_t ctx)
{
    fmpq_mpoly_one(res, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpq_mpoly_gens(gr_vec_t res, gr_ctx_t ctx)
{
    slong i, n;

    n = MPOLYNOMIAL_MCTX(ctx)->zctx->minfo->nvars;

    gr_vec_set_length(res, n, ctx);
    for (i = 0; i < n; i++)
        fmpq_mpoly_gen(((fmpq_mpoly_struct *) res->entries) + i, i, MPOLYNOMIAL_MCTX(ctx));

    return GR_SUCCESS;
}

static int
_gr_fmpq_mpoly_set(fmpq_mpoly_t res, const fmpq_mpoly_t mat, gr_ctx_t ctx)
{
    fmpq_mpoly_set(res, mat, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpq_mpoly_set_si(fmpq_mpoly_t res, slong v, gr_ctx_t ctx)
{
    fmpq_mpoly_set_si(res, v, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpq_mpoly_set_ui(fmpq_mpoly_t res, ulong v, gr_ctx_t ctx)
{
    fmpq_mpoly_set_ui(res, v, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpq_mpoly_set_fmpz(fmpq_mpoly_t res, const fmpz_t v, gr_ctx_t ctx)
{
    fmpq_mpoly_set_fmpz(res, v, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpq_mpoly_set_fmpq(fmpq_mpoly_t res, const fmpq_t v, gr_ctx_t ctx)
{
    fmpq_mpoly_set_fmpq(res, v, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpq_mpoly_set_other(fmpq_mpoly_t res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
{
    if (x_ctx->which_ring == GR_CTX_FMPQ_MPOLY)
    {
        if (MPOLYNOMIAL_MCTX(ctx)->zctx->minfo->nvars == MPOLYNOMIAL_MCTX(x_ctx)->zctx->minfo->nvars &&
            MPOLYNOMIAL_MCTX(ctx)->zctx->minfo->ord == MPOLYNOMIAL_MCTX(x_ctx)->zctx->minfo->ord)
        {
            fmpq_mpoly_set(res, x, MPOLYNOMIAL_MCTX(ctx));
            return GR_SUCCESS;
        }
    }

    return gr_generic_set_other(res, x, x_ctx, ctx);
}

static int
_gr_fmpq_mpoly_neg(fmpq_mpoly_t res, const fmpq_mpoly_t mat, gr_ctx_t ctx)
{
    fmpq_mpoly_neg(res, mat, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpq_mpoly_add(fmpq_mpoly_t res, const fmpq_mpoly_t poly1, const fmpq_mpoly_t poly2, gr_ctx_t ctx)
{
    if (poly1->zpoly->length + poly2->zpoly->length > ctx->size_limit)
    {
        fmpq_mpoly_zero(res, MPOLYNOMIAL_MCTX(ctx));
        return GR_UNABLE;
    }

    fmpq_mpoly_add(res, poly1, poly2, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpq_mpoly_add_si(fmpq_mpoly_t res, const fmpq_mpoly_t poly1, slong c, gr_ctx_t ctx)
{
    fmpq_mpoly_add_si(res, poly1, c, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpq_mpoly_add_ui(fmpq_mpoly_t res, const fmpq_mpoly_t poly1, ulong c, gr_ctx_t ctx)
{
    fmpq_mpoly_add_ui(res, poly1, c, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpq_mpoly_add_fmpz(fmpq_mpoly_t res, const fmpq_mpoly_t poly1, const fmpz_t c, gr_ctx_t ctx)
{
    fmpq_mpoly_add_fmpz(res, poly1, c, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpq_mpoly_sub(fmpq_mpoly_t res, const fmpq_mpoly_t poly1, const fmpq_mpoly_t poly2, gr_ctx_t ctx)
{
    if (poly1->zpoly->length + poly2->zpoly->length > ctx->size_limit)
    {
        fmpq_mpoly_zero(res, MPOLYNOMIAL_MCTX(ctx));
        return GR_UNABLE;
    }

    fmpq_mpoly_sub(res, poly1, poly2, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpq_mpoly_sub_si(fmpq_mpoly_t res, const fmpq_mpoly_t poly1, slong c, gr_ctx_t ctx)
{
    fmpq_mpoly_sub_si(res, poly1, c, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpq_mpoly_sub_ui(fmpq_mpoly_t res, const fmpq_mpoly_t poly1, ulong c, gr_ctx_t ctx)
{
    fmpq_mpoly_sub_ui(res, poly1, c, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpq_mpoly_sub_fmpz(fmpq_mpoly_t res, const fmpq_mpoly_t poly1, const fmpz_t c, gr_ctx_t ctx)
{
    fmpq_mpoly_sub_fmpz(res, poly1, c, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpq_mpoly_mul(fmpq_mpoly_t res, const fmpq_mpoly_t poly1, const fmpq_mpoly_t poly2, gr_ctx_t ctx)
{
    if (poly1->zpoly->length * poly2->zpoly->length > ctx->size_limit) /* todo: * can overflow */
    {
        fmpq_mpoly_zero(res, MPOLYNOMIAL_MCTX(ctx));
        return GR_UNABLE;
    }

    fmpq_mpoly_mul(res, poly1, poly2, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpq_mpoly_mul_si(fmpq_mpoly_t res, const fmpq_mpoly_t poly1, slong c, gr_ctx_t ctx)
{
    fmpq_mpoly_scalar_mul_si(res, poly1, c, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpq_mpoly_mul_ui(fmpq_mpoly_t res, const fmpq_mpoly_t poly1, ulong c, gr_ctx_t ctx)
{
    fmpq_mpoly_scalar_mul_ui(res, poly1, c, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpq_mpoly_mul_fmpz(fmpq_mpoly_t res, const fmpq_mpoly_t poly1, const fmpz_t c, gr_ctx_t ctx)
{
    fmpq_mpoly_scalar_mul_fmpz(res, poly1, c, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpq_mpoly_div(fmpq_mpoly_t res, const fmpq_mpoly_t poly1, const fmpq_mpoly_t poly2, gr_ctx_t ctx)
{
    if (poly2->zpoly->length == 0)
        return GR_DOMAIN;

    if (fmpq_mpoly_divides(res, poly1, poly2, MPOLYNOMIAL_MCTX(ctx)))
        return GR_SUCCESS;
    else
        return GR_DOMAIN;
}

static int
_gr_fmpq_mpoly_divexact(fmpq_mpoly_t res, const fmpq_mpoly_t poly1, const fmpq_mpoly_t poly2, gr_ctx_t ctx)
{
    if (poly2->zpoly->length == 0)
        return GR_DOMAIN;

    fmpq_mpoly_div(res, poly1, poly2, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpq_mpoly_divexact_si(fmpq_mpoly_t res, const fmpq_mpoly_t poly1, slong c, gr_ctx_t ctx)
{
    if (c == 0)
        return GR_DOMAIN;

    fmpq_mpoly_scalar_div_si(res, poly1, c, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpq_mpoly_divexact_ui(fmpq_mpoly_t res, const fmpq_mpoly_t poly1, ulong c, gr_ctx_t ctx)
{
    if (c == 0)
        return GR_DOMAIN;

    fmpq_mpoly_scalar_div_ui(res, poly1, c, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpq_mpoly_divexact_fmpz(fmpq_mpoly_t res, const fmpq_mpoly_t poly1, const fmpz_t c, gr_ctx_t ctx)
{
    if (fmpz_is_zero(c))
        return GR_DOMAIN;

    fmpq_mpoly_scalar_div_fmpz(res, poly1, c, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static truth_t
_gr_fmpq_mpoly_is_invertible(const fmpq_mpoly_t c, gr_ctx_t ctx)
{
    if (fmpq_mpoly_is_fmpq(c, MPOLYNOMIAL_MCTX(ctx)) && !fmpq_mpoly_is_zero(c, MPOLYNOMIAL_MCTX(ctx)))
        return T_TRUE;

    return T_FALSE;
}

static int
_gr_fmpq_mpoly_inv(fmpq_mpoly_t res, const fmpq_mpoly_t c, gr_ctx_t ctx)
{
    if (fmpq_mpoly_is_fmpq(c, MPOLYNOMIAL_MCTX(ctx)) && !fmpq_mpoly_is_zero(c, MPOLYNOMIAL_MCTX(ctx)))
    {
        /* todo: optimize */
        fmpq_t t;
        fmpq_init(t);
        fmpq_mpoly_get_fmpq(t, c, MPOLYNOMIAL_MCTX(ctx));
        fmpq_inv(t, t);
        fmpq_mpoly_set_fmpq(res, t, MPOLYNOMIAL_MCTX(ctx));
        fmpq_clear(t);
        return GR_SUCCESS;
    }

    return GR_DOMAIN;
}

static int
_gr_fmpq_mpoly_pow_ui(fmpq_mpoly_t res, const fmpq_mpoly_t poly1, ulong c, gr_ctx_t ctx)
{
    /* todo: size limit */

    if (fmpq_mpoly_pow_ui(res, poly1, c, MPOLYNOMIAL_MCTX(ctx)))
        return GR_SUCCESS;
    else
        return GR_UNABLE;
}

static int
_gr_fmpq_mpoly_pow_fmpz(fmpq_mpoly_t res, const fmpq_mpoly_t poly1, const fmpz_t c, gr_ctx_t ctx)
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
            status = _gr_fmpq_mpoly_pow_fmpz(res, res, e, ctx);
            fmpz_clear(e);
        }

        return status;
    }

    if (fmpq_mpoly_pow_fmpz(res, poly1, c, MPOLYNOMIAL_MCTX(ctx)))
        return GR_SUCCESS;
    else
        return GR_UNABLE;
}

static int
_gr_fmpq_mpoly_derivative_gen(fmpq_mpoly_t res, const fmpq_mpoly_t poly, slong var, gr_ctx_t ctx)
{
    if (var >= 0 && var < MPOLYNOMIAL_MCTX(ctx)->zctx->minfo->nvars)
    {
        fmpq_mpoly_derivative(res, poly, var, MPOLYNOMIAL_MCTX(ctx));
        return GR_SUCCESS;
    }
    return GR_DOMAIN;
}

static truth_t
_gr_fmpq_mpoly_is_square(const fmpq_mpoly_t poly, gr_ctx_t ctx)
{
    return fmpq_mpoly_is_square(poly, MPOLYNOMIAL_MCTX(ctx)) ? T_TRUE : T_FALSE;
}

static int
_gr_fmpq_mpoly_sqrt(fmpq_mpoly_t res, const fmpq_mpoly_t poly, gr_ctx_t ctx)
{
    if (fmpq_mpoly_sqrt(res, poly, MPOLYNOMIAL_MCTX(ctx)))
        return GR_SUCCESS;
    else
        return GR_DOMAIN;
}

static int
_gr_fmpq_mpoly_gcd(fmpq_mpoly_t res, const fmpq_mpoly_t poly1, const fmpq_mpoly_t poly2, gr_ctx_t ctx)
{
    if (fmpq_mpoly_gcd(res, poly1, poly2, MPOLYNOMIAL_MCTX(ctx)))
        return GR_SUCCESS;

    return GR_UNABLE;
}

static int
_gr_fmpq_mpoly_denominator(fmpq_mpoly_t res, const fmpq_mpoly_t poly, gr_ctx_t ctx)
{
    fmpq_set_fmpz(res->content, fmpq_denref(poly->content));
    fmpz_mpoly_one(res->zpoly, MPOLYNOMIAL_MCTX(ctx)->zctx);
    return GR_SUCCESS;
}

static int
_gr_fmpq_mpoly_numerator(fmpq_mpoly_t res, const fmpq_mpoly_t poly, gr_ctx_t ctx)
{
    fmpq_mpoly_set(res, poly, MPOLYNOMIAL_MCTX(ctx));
    fmpz_one(fmpq_denref(res->content));
    return GR_SUCCESS;
}

static int
_gr_fmpq_mpoly_factor(fmpq_mpoly_t c, gr_vec_t factors, gr_vec_t exponents, gr_srcptr x, int flags, gr_ctx_t ctx)
{
    fmpq_mpoly_factor_t fac;
    gr_ctx_t ZZ;
    slong i;
    int status = GR_SUCCESS;

    fmpq_mpoly_factor_init(fac, MPOLYNOMIAL_MCTX(ctx));

    if (fmpq_mpoly_factor(fac, x, MPOLYNOMIAL_MCTX(ctx)))
    {
        fmpq_mpoly_set_fmpq(c, fac->constant, MPOLYNOMIAL_MCTX(ctx));

        gr_ctx_init_fmpz(ZZ);

        gr_vec_set_length(factors, fac->num, ctx);
        gr_vec_set_length(exponents, fac->num, ZZ);

        for (i = 0; i < fac->num; i++)
        {
            fmpq_mpoly_swap((fmpq_mpoly_struct *) (factors->entries) + i, fac->poly + i, MPOLYNOMIAL_MCTX(ctx));
            fmpz_swap((fmpz *) (exponents->entries) + i, fac->exp + i);
        }

        gr_ctx_clear(ZZ);
    }
    else
    {
        status = GR_UNABLE;
    }

    fmpq_mpoly_factor_clear(fac, MPOLYNOMIAL_MCTX(ctx));

    return status;
}


int _gr_fmpq_mpoly_methods_initialized = 0;

gr_static_method_table _gr_fmpq_mpoly_methods;

gr_method_tab_input _gr_fmpq_mpoly_methods_input[] =
{
    {GR_METHOD_CTX_WRITE,   (gr_funcptr) _gr_fmpq_mpoly_ctx_write},
    {GR_METHOD_CTX_CLEAR,   (gr_funcptr) _gr_fmpq_mpoly_ctx_clear},
    {GR_METHOD_CTX_IS_RING,                     (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_COMMUTATIVE_RING,         (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_INTEGRAL_DOMAIN,          (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_UNIQUE_FACTORIZATION_DOMAIN,    (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_FIELD,                    (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_FINITE,                   (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_FINITE_CHARACTERISTIC,    (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_RATIONAL_VECTOR_SPACE, (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_REAL_VECTOR_SPACE, (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_COMPLEX_VECTOR_SPACE, (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_THREADSAFE,               (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_SET_GEN_NAMES,               (gr_funcptr) _gr_fmpq_mpoly_ctx_set_gen_names},
    {GR_METHOD_CTX_NGENS,   (gr_funcptr) _gr_fmpq_mpoly_ctx_ngens},
    {GR_METHOD_CTX_GEN_NAME, (gr_funcptr) _gr_fmpq_mpoly_ctx_gen_name},
    {GR_METHOD_INIT,        (gr_funcptr) _gr_fmpq_mpoly_init},
    {GR_METHOD_CLEAR,       (gr_funcptr) _gr_fmpq_mpoly_clear},
    {GR_METHOD_SWAP,        (gr_funcptr) _gr_fmpq_mpoly_swap},
    {GR_METHOD_SET_SHALLOW, (gr_funcptr) _gr_fmpq_mpoly_set_shallow},
    {GR_METHOD_RANDTEST,    (gr_funcptr) _gr_fmpq_mpoly_randtest},
    {GR_METHOD_RANDTEST_SMALL,    (gr_funcptr) _gr_fmpq_mpoly_randtest_small},
    {_GR_METHOD_LENGTH,     (gr_funcptr) _gr_fmpq_mpoly_length},
    {GR_METHOD_WRITE,       (gr_funcptr) _gr_fmpq_mpoly_write},
    {GR_METHOD_ZERO,        (gr_funcptr) _gr_fmpq_mpoly_zero},
    {GR_METHOD_ONE,         (gr_funcptr) _gr_fmpq_mpoly_one},
    {GR_METHOD_IS_ZERO,     (gr_funcptr) _gr_fmpq_mpoly_is_zero},
    {GR_METHOD_IS_ONE,      (gr_funcptr) _gr_fmpq_mpoly_is_one},
    {GR_METHOD_GENS,        (gr_funcptr) _gr_fmpq_mpoly_gens},
    {GR_METHOD_EQUAL,       (gr_funcptr) _gr_fmpq_mpoly_equal},
    {GR_METHOD_SET,         (gr_funcptr) _gr_fmpq_mpoly_set},
    {GR_METHOD_SET_UI,      (gr_funcptr) _gr_fmpq_mpoly_set_ui},
    {GR_METHOD_SET_SI,      (gr_funcptr) _gr_fmpq_mpoly_set_si},
    {GR_METHOD_SET_FMPZ,    (gr_funcptr) _gr_fmpq_mpoly_set_fmpz},
    {GR_METHOD_SET_FMPQ,    (gr_funcptr) _gr_fmpq_mpoly_set_fmpq},
    {GR_METHOD_SET_OTHER,   (gr_funcptr) _gr_fmpq_mpoly_set_other},
    {GR_METHOD_SET_STR,     (gr_funcptr) gr_generic_set_str_balance_additions},
    {GR_METHOD_NEG,         (gr_funcptr) _gr_fmpq_mpoly_neg},
    {GR_METHOD_ADD,         (gr_funcptr) _gr_fmpq_mpoly_add},
    {GR_METHOD_ADD_SI,      (gr_funcptr) _gr_fmpq_mpoly_add_si},
    {GR_METHOD_ADD_UI,      (gr_funcptr) _gr_fmpq_mpoly_add_ui},
    {GR_METHOD_ADD_FMPZ,    (gr_funcptr) _gr_fmpq_mpoly_add_fmpz},
    {GR_METHOD_SUB,         (gr_funcptr) _gr_fmpq_mpoly_sub},
    {GR_METHOD_SUB_SI,      (gr_funcptr) _gr_fmpq_mpoly_sub_si},
    {GR_METHOD_SUB_UI,      (gr_funcptr) _gr_fmpq_mpoly_sub_ui},
    {GR_METHOD_SUB_FMPZ,    (gr_funcptr) _gr_fmpq_mpoly_sub_fmpz},
    {GR_METHOD_MUL,         (gr_funcptr) _gr_fmpq_mpoly_mul},
    {GR_METHOD_MUL_SI,      (gr_funcptr) _gr_fmpq_mpoly_mul_si},
    {GR_METHOD_MUL_UI,      (gr_funcptr) _gr_fmpq_mpoly_mul_ui},
    {GR_METHOD_MUL_FMPZ,    (gr_funcptr) _gr_fmpq_mpoly_mul_fmpz},
    {GR_METHOD_DIV,         (gr_funcptr) _gr_fmpq_mpoly_div},
    {GR_METHOD_DIVEXACT,         (gr_funcptr) _gr_fmpq_mpoly_divexact},
    {GR_METHOD_DIVEXACT_SI,      (gr_funcptr) _gr_fmpq_mpoly_divexact_si},
    {GR_METHOD_DIVEXACT_UI,      (gr_funcptr) _gr_fmpq_mpoly_divexact_ui},
    {GR_METHOD_DIVEXACT_FMPZ,    (gr_funcptr) _gr_fmpq_mpoly_divexact_fmpz},
    {GR_METHOD_INV,             (gr_funcptr) _gr_fmpq_mpoly_inv},
    {GR_METHOD_IS_INVERTIBLE,   (gr_funcptr) _gr_fmpq_mpoly_is_invertible},
    {GR_METHOD_POW_UI,      (gr_funcptr) _gr_fmpq_mpoly_pow_ui},
    {GR_METHOD_POW_FMPZ,    (gr_funcptr) _gr_fmpq_mpoly_pow_fmpz},
    {GR_METHOD_DERIVATIVE_GEN,  (gr_funcptr) _gr_fmpq_mpoly_derivative_gen},
    {GR_METHOD_SQRT,        (gr_funcptr) _gr_fmpq_mpoly_sqrt},
    {GR_METHOD_IS_SQUARE,   (gr_funcptr) _gr_fmpq_mpoly_is_square},
    {GR_METHOD_NUMERATOR,    (gr_funcptr) _gr_fmpq_mpoly_numerator},
    {GR_METHOD_DENOMINATOR,  (gr_funcptr) _gr_fmpq_mpoly_denominator},
    {GR_METHOD_GCD,         (gr_funcptr) _gr_fmpq_mpoly_gcd},
    {GR_METHOD_FACTOR,      (gr_funcptr) _gr_fmpq_mpoly_factor},
    {0,                     (gr_funcptr) NULL},
};

void
gr_ctx_init_fmpq_mpoly(gr_ctx_t ctx, slong nvars, const ordering_t ord)
{
    ctx->which_ring = GR_CTX_FMPQ_MPOLY;
    ctx->sizeof_elem = sizeof(fmpq_mpoly_struct);
    GR_CTX_DATA_AS_PTR(ctx) = flint_malloc(sizeof(_gr_fmpq_mpoly_ctx_t));
    ctx->size_limit = WORD_MAX;

    fmpq_mpoly_ctx_init(MPOLYNOMIAL_MCTX(ctx), nvars, ord);
    MPOLYNOMIAL_CTX(ctx)->vars = NULL;

    ctx->methods = _gr_fmpq_mpoly_methods;

    if (!_gr_fmpq_mpoly_methods_initialized)
    {
        gr_method_tab_init(_gr_fmpq_mpoly_methods, _gr_fmpq_mpoly_methods_input);
        _gr_fmpq_mpoly_methods_initialized = 1;
    }
}
