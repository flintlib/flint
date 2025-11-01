/*
    Copyright (C) 2023 Fredrik Johansson
    Copyright (C) 2025 Andrii Yanovets

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "fmpq.h"
#include "fmpz_mod_mpoly_q.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_generic.h"

typedef struct
{
    fmpz_mod_mpoly_ctx_t mctx;
    char ** vars;
}
_gr_fmpz_mod_mpoly_ctx_t;

#define MPOLYNOMIAL_CTX(ring_ctx) ((_gr_fmpz_mod_mpoly_ctx_t *)(GR_CTX_DATA_AS_PTR(ring_ctx)))
#define MPOLYNOMIAL_MCTX(ring_ctx) (MPOLYNOMIAL_CTX(ring_ctx)->mctx)

static int _gr_fmpz_mod_mpoly_q_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    gr_stream_write(out, "Fraction field of multivariate polynomials over finite field (fmpz_mod) mod ");
    gr_stream_write_fmpz(out, MPOLYNOMIAL_MCTX(ctx)->ffinfo->n);
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

static void _gr_fmpz_mod_mpoly_q_ctx_clear(gr_ctx_t ctx)
{
    if (MPOLYNOMIAL_CTX(ctx)->vars != NULL)
    {
        slong i;
        for (i = 0; i < MPOLYNOMIAL_MCTX(ctx)->minfo->nvars; i++)
            flint_free(MPOLYNOMIAL_CTX(ctx)->vars[i]);
        flint_free(MPOLYNOMIAL_CTX(ctx)->vars);
    }

    fmpz_mod_mpoly_ctx_clear(MPOLYNOMIAL_MCTX(ctx));
    flint_free(GR_CTX_DATA_AS_PTR(ctx));  
}

static int
_gr_fmpz_mod_mpoly_q_ctx_set_gen_names(gr_ctx_t ctx, const char ** s)
{
    slong i, nvars, len;

    nvars = MPOLYNOMIAL_MCTX(ctx)->minfo->nvars;

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
_gr_fmpz_mod_mpoly_q_ctx_ngens(slong * ngens, gr_ctx_t ctx)
{
    * ngens = MPOLYNOMIAL_MCTX(ctx)->minfo->nvars;
    return GR_SUCCESS;
}

static int
_gr_fmpz_mod_mpoly_q_ctx_gen_name(char ** name, slong i, gr_ctx_t ctx)
{
    if (i < 0 || i >= MPOLYNOMIAL_MCTX(ctx)->minfo->nvars)
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
_gr_fmpz_mod_mpoly_q_init(fmpz_mod_mpoly_q_t res, gr_ctx_t ctx)
{
    fmpz_mod_mpoly_q_init(res, MPOLYNOMIAL_MCTX(ctx));
}

static void
_gr_fmpz_mod_mpoly_q_clear(fmpz_mod_mpoly_q_t res, gr_ctx_t ctx)
{
    fmpz_mod_mpoly_q_clear(res, MPOLYNOMIAL_MCTX(ctx));
}

static void
_gr_fmpz_mod_mpoly_q_swap(fmpz_mod_mpoly_q_t poly1, fmpz_mod_mpoly_q_t poly2, gr_ctx_t ctx)
{
    fmpz_mod_mpoly_q_swap(poly1, poly2, MPOLYNOMIAL_MCTX(ctx));
}

static void
_gr_fmpz_mod_mpoly_q_set_shallow(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t poly, gr_ctx_t ctx)
{
    *res = *poly;
}

static int
_gr_fmpz_mod_mpoly_q_randtest(fmpz_mod_mpoly_q_t res, flint_rand_t state, gr_ctx_t ctx)
{
    fmpz_mod_mpoly_q_randtest(res, state, n_randint(state, 5), 1 + n_randint(state, 3), MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpz_mod_mpoly_q_randtest_small(fmpz_mod_mpoly_q_t res, flint_rand_t state, gr_ctx_t ctx)
{
    fmpz_mod_mpoly_q_randtest(res, state, n_randint(state, 3), 1 + n_randint(state, 3), MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static slong
_gr_fmpz_mod_mpoly_q_length(const fmpz_mod_mpoly_q_t x, gr_ctx_t ctx)
{
    return fmpz_mod_mpoly_q_numref(x)->length + fmpz_mod_mpoly_q_denref(x)->length;
}

static int
_gr_fmpz_mod_mpoly_q_write(gr_stream_t out, fmpz_mod_mpoly_q_t f, gr_ctx_t ctx)
{
    if (fmpz_mod_mpoly_is_one(fmpz_mod_mpoly_q_denref(f), MPOLYNOMIAL_MCTX(ctx)))
    {
        gr_stream_write_free(out, fmpz_mod_mpoly_get_str_pretty(fmpz_mod_mpoly_q_numref(f), (const char **) MPOLYNOMIAL_CTX(ctx)->vars, MPOLYNOMIAL_MCTX(ctx)));
    }
    else if (fmpz_mod_mpoly_is_fmpz(fmpz_mod_mpoly_q_denref(f), MPOLYNOMIAL_MCTX(ctx)))
    {
        gr_stream_write(out, "(");
        gr_stream_write_free(out, fmpz_mod_mpoly_get_str_pretty(fmpz_mod_mpoly_q_numref(f), (const char **) MPOLYNOMIAL_CTX(ctx)->vars, MPOLYNOMIAL_MCTX(ctx)));
        gr_stream_write(out, ")/");
        gr_stream_write_free(out, fmpz_mod_mpoly_get_str_pretty(fmpz_mod_mpoly_q_denref(f), (const char **) MPOLYNOMIAL_CTX(ctx)->vars, MPOLYNOMIAL_MCTX(ctx)));
    }
    else
    {
        gr_stream_write(out, "(");
        gr_stream_write_free(out, fmpz_mod_mpoly_get_str_pretty(fmpz_mod_mpoly_q_numref(f), (const char **) MPOLYNOMIAL_CTX(ctx)->vars, MPOLYNOMIAL_MCTX(ctx)));
        gr_stream_write(out, ")/(");
        gr_stream_write_free(out, fmpz_mod_mpoly_get_str_pretty(fmpz_mod_mpoly_q_denref(f), (const char **) MPOLYNOMIAL_CTX(ctx)->vars, MPOLYNOMIAL_MCTX(ctx)));
        gr_stream_write(out, ")");
    }

    return GR_SUCCESS;
}

static truth_t
_gr_fmpz_mod_mpoly_q_equal(const fmpz_mod_mpoly_q_t poly1, const fmpz_mod_mpoly_q_t poly2, gr_ctx_t ctx)
{
    return fmpz_mod_mpoly_q_equal(poly1, poly2, MPOLYNOMIAL_MCTX(ctx)) ? T_TRUE : T_FALSE;
}

static truth_t
_gr_fmpz_mod_mpoly_q_is_zero(const fmpz_mod_mpoly_q_t poly, gr_ctx_t ctx)
{
    return fmpz_mod_mpoly_q_is_zero(poly, MPOLYNOMIAL_MCTX(ctx)) ? T_TRUE : T_FALSE;
}

static truth_t
_gr_fmpz_mod_mpoly_q_is_one(const fmpz_mod_mpoly_q_t poly, gr_ctx_t ctx)
{
    return fmpz_mod_mpoly_q_is_one(poly, MPOLYNOMIAL_MCTX(ctx)) ? T_TRUE : T_FALSE;
}

static int
_gr_fmpz_mod_mpoly_q_zero(fmpz_mod_mpoly_q_t res, gr_ctx_t ctx)
{
    fmpz_mod_mpoly_q_zero(res, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpz_mod_mpoly_q_one(fmpz_mod_mpoly_q_t res, gr_ctx_t ctx)
{
    fmpz_mod_mpoly_q_one(res, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpz_mod_mpoly_q_gens(gr_vec_t res, gr_ctx_t ctx)
{
    slong i, n;

    n = MPOLYNOMIAL_MCTX(ctx)->minfo->nvars;

    gr_vec_set_length(res, n, ctx);
    for (i = 0; i < n; i++)
        fmpz_mod_mpoly_q_gen(((fmpz_mod_mpoly_q_struct *) res->entries) + i, i, MPOLYNOMIAL_MCTX(ctx));

    return GR_SUCCESS;
}

static int
_gr_fmpz_mod_mpoly_q_set(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t mat, gr_ctx_t ctx)
{
    fmpz_mod_mpoly_q_set(res, mat, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpz_mod_mpoly_q_set_si(fmpz_mod_mpoly_q_t res, slong v, gr_ctx_t ctx)
{
    fmpz_mod_mpoly_q_set_si(res, v, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpz_mod_mpoly_q_set_ui(fmpz_mod_mpoly_q_t res, ulong v, gr_ctx_t ctx)
{
    fmpz_mod_mpoly_set_ui(fmpz_mod_mpoly_q_numref(res), v, MPOLYNOMIAL_MCTX(ctx));
    fmpz_mod_mpoly_one(fmpz_mod_mpoly_q_denref(res), MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpz_mod_mpoly_q_set_fmpz(fmpz_mod_mpoly_q_t res, const fmpz_t v, gr_ctx_t ctx)
{
    fmpz_mod_mpoly_q_set_fmpz(res, v, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpz_mod_mpoly_q_set_fmpq(fmpz_mod_mpoly_q_t res, const fmpq_t v, gr_ctx_t ctx)
{
    return fmpz_mod_mpoly_q_set_fmpq(res, v, MPOLYNOMIAL_MCTX(ctx)) ? GR_SUCCESS : GR_DOMAIN;
}

static int
_gr_fmpz_mod_mpoly_q_neg(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t mat, gr_ctx_t ctx)
{
    fmpz_mod_mpoly_q_neg(res, mat, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpz_mod_mpoly_q_add(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t poly1, const fmpz_mod_mpoly_q_t poly2, gr_ctx_t ctx)
{
    fmpz_mod_mpoly_q_add(res, poly1, poly2, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpz_mod_mpoly_q_add_si(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t poly1, slong c, gr_ctx_t ctx)
{
    fmpz_mod_mpoly_q_add_si(res, poly1, c, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}



static int
_gr_fmpz_mod_mpoly_q_add_fmpz(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t poly1, const fmpz_t c, gr_ctx_t ctx)
{
    fmpz_mod_mpoly_q_add_fmpz(res, poly1, c, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpz_mod_mpoly_q_add_fmpq(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t poly1, const fmpq_t c, gr_ctx_t ctx)
{
    return fmpz_mod_mpoly_q_add_fmpq(res, poly1, c, MPOLYNOMIAL_MCTX(ctx)) ? GR_SUCCESS : GR_DOMAIN;
}

static int
_gr_fmpz_mod_mpoly_q_sub(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t poly1, const fmpz_mod_mpoly_q_t poly2, gr_ctx_t ctx)
{
    fmpz_mod_mpoly_q_sub(res, poly1, poly2, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpz_mod_mpoly_q_sub_si(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t poly1, slong c, gr_ctx_t ctx)
{
    fmpz_mod_mpoly_q_sub_si(res, poly1, c, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpz_mod_mpoly_q_sub_fmpz(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t poly1, const fmpz_t c, gr_ctx_t ctx)
{
    fmpz_mod_mpoly_q_sub_fmpz(res, poly1, c, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpz_mod_mpoly_q_sub_fmpq(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t poly1, const fmpq_t c, gr_ctx_t ctx)
{
    return fmpz_mod_mpoly_q_sub_fmpq(res, poly1, c, MPOLYNOMIAL_MCTX(ctx)) ? GR_SUCCESS : GR_DOMAIN;
}

static int
_gr_fmpz_mod_mpoly_q_mul(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t poly1, const fmpz_mod_mpoly_q_t poly2, gr_ctx_t ctx)
{
    fmpz_mod_mpoly_q_mul(res, poly1, poly2, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpz_mod_mpoly_q_mul_si(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t poly1, slong c, gr_ctx_t ctx)
{
    fmpz_mod_mpoly_q_mul_si(res, poly1, c, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpz_mod_mpoly_q_mul_fmpz(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t poly1, const fmpz_t c, gr_ctx_t ctx)
{
    fmpz_mod_mpoly_q_mul_fmpz(res, poly1, c, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpz_mod_mpoly_q_mul_fmpq(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t poly1, const fmpq_t c, gr_ctx_t ctx)
{
    return fmpz_mod_mpoly_q_mul_fmpq(res, poly1, c, MPOLYNOMIAL_MCTX(ctx)) ? GR_SUCCESS : GR_DOMAIN;
}

static int
_gr_fmpz_mod_mpoly_q_div(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t poly1, const fmpz_mod_mpoly_q_t poly2, gr_ctx_t ctx)
{
    if (fmpz_mod_mpoly_q_is_zero(poly2, MPOLYNOMIAL_MCTX(ctx)))
        return GR_DOMAIN;

    fmpz_mod_mpoly_q_div(res, poly1, poly2, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpz_mod_mpoly_q_div_si(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t poly1, slong c, gr_ctx_t ctx)
{
    return fmpz_mod_mpoly_q_div_si(res, poly1, c, MPOLYNOMIAL_MCTX(ctx)) ? GR_SUCCESS : GR_DOMAIN;
}

static int
_gr_fmpz_mod_mpoly_q_div_fmpz(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t poly1, const fmpz_t c, gr_ctx_t ctx)
{
    return fmpz_mod_mpoly_q_div_fmpz(res, poly1, c, MPOLYNOMIAL_MCTX(ctx)) ? GR_SUCCESS : GR_DOMAIN;
}

static int
_gr_fmpz_mod_mpoly_q_div_fmpq(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t poly1, const fmpq_t c, gr_ctx_t ctx)
{
    return fmpz_mod_mpoly_q_div_fmpq(res, poly1, c, MPOLYNOMIAL_MCTX(ctx)) ? GR_SUCCESS : GR_DOMAIN;
}

static truth_t
_gr_fmpz_mod_mpoly_q_is_invertible(const fmpz_mod_mpoly_q_t c, gr_ctx_t ctx)
{
    return fmpz_mod_mpoly_q_is_zero(c, MPOLYNOMIAL_MCTX(ctx)) ? T_FALSE : T_TRUE;
}

static int
_gr_fmpz_mod_mpoly_q_inv(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t c, gr_ctx_t ctx)
{
    if (!fmpz_mod_mpoly_q_is_zero(c, MPOLYNOMIAL_MCTX(ctx)))
    {
        fmpz_mod_mpoly_q_inv(res, c, MPOLYNOMIAL_MCTX(ctx));
        return GR_SUCCESS;
    }

    return GR_DOMAIN;
}

static int
_gr_fmpz_mod_mpoly_q_pow_ui(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t poly1, ulong c, gr_ctx_t ctx)
{
    if (fmpz_mod_mpoly_pow_ui(fmpz_mod_mpoly_q_numref(res), fmpz_mod_mpoly_q_numref(poly1), c, MPOLYNOMIAL_MCTX(ctx)) &&
        fmpz_mod_mpoly_pow_ui(fmpz_mod_mpoly_q_denref(res), fmpz_mod_mpoly_q_denref(poly1), c, MPOLYNOMIAL_MCTX(ctx)))
        return GR_SUCCESS;
    else
        return GR_UNABLE;
}

static int
_gr_fmpz_mod_mpoly_q_pow_fmpz(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t poly1, const fmpz_t c, gr_ctx_t ctx)
{
    if (fmpz_sgn(c) < 0)
    {
        int status;

        status = gr_inv(res, poly1, ctx);

        if (status == GR_SUCCESS)
        {
            fmpz_t e;
            fmpz_init(e);
            fmpz_neg(e, c);
            status = _gr_fmpz_mod_mpoly_q_pow_fmpz(res, res, e, ctx);
            fmpz_clear(e);
        }

        return status;
    }

    if (fmpz_mod_mpoly_pow_fmpz(fmpz_mod_mpoly_q_numref(res), fmpz_mod_mpoly_q_numref(poly1), c, MPOLYNOMIAL_MCTX(ctx)) &&
        fmpz_mod_mpoly_pow_fmpz(fmpz_mod_mpoly_q_denref(res), fmpz_mod_mpoly_q_denref(poly1), c, MPOLYNOMIAL_MCTX(ctx)))
        return GR_SUCCESS;
    else
        return GR_UNABLE;

}

static int
_gr_fmpz_mod_mpoly_q_numerator(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const gr_ctx_t ctx)
{
    fmpz_mod_mpoly_set(fmpz_mod_mpoly_q_numref(res), fmpz_mod_mpoly_q_numref(x), MPOLYNOMIAL_MCTX(ctx));
    fmpz_mod_mpoly_one(fmpz_mod_mpoly_q_denref(res), MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpz_mod_mpoly_q_denominator(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const gr_ctx_t ctx)
{
    fmpz_mod_mpoly_set(fmpz_mod_mpoly_q_numref(res), fmpz_mod_mpoly_q_denref(x), MPOLYNOMIAL_MCTX(ctx));
    fmpz_mod_mpoly_one(fmpz_mod_mpoly_q_denref(res), MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

int _gr_fmpz_mod_mpoly_q_methods_initialized = 0;

gr_static_method_table _gr_fmpz_mod_mpoly_q_methods;

gr_method_tab_input _gr_fmpz_mod_mpoly_q_methods_input[] =
{
    {GR_METHOD_CTX_WRITE,   (gr_funcptr) _gr_fmpz_mod_mpoly_q_ctx_write},
    {GR_METHOD_CTX_CLEAR,   (gr_funcptr) _gr_fmpz_mod_mpoly_q_ctx_clear},
    {GR_METHOD_CTX_IS_RING,                     (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_COMMUTATIVE_RING,         (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_INTEGRAL_DOMAIN,          (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_UNIQUE_FACTORIZATION_DOMAIN,    (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_FIELD,                    (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_FINITE,                   (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_FINITE_CHARACTERISTIC,    (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_THREADSAFE,               (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_SET_GEN_NAMES,               (gr_funcptr) _gr_fmpz_mod_mpoly_q_ctx_set_gen_names},
    {GR_METHOD_CTX_NGENS,   (gr_funcptr) _gr_fmpz_mod_mpoly_q_ctx_ngens},
    {GR_METHOD_CTX_GEN_NAME, (gr_funcptr) _gr_fmpz_mod_mpoly_q_ctx_gen_name},
    {GR_METHOD_INIT,        (gr_funcptr) _gr_fmpz_mod_mpoly_q_init},
    {GR_METHOD_CLEAR,       (gr_funcptr) _gr_fmpz_mod_mpoly_q_clear},
    {GR_METHOD_SWAP,        (gr_funcptr) _gr_fmpz_mod_mpoly_q_swap},
    {GR_METHOD_SET_SHALLOW, (gr_funcptr) _gr_fmpz_mod_mpoly_q_set_shallow},
    {GR_METHOD_RANDTEST,    (gr_funcptr) _gr_fmpz_mod_mpoly_q_randtest},
    {GR_METHOD_RANDTEST_SMALL,    (gr_funcptr) _gr_fmpz_mod_mpoly_q_randtest_small},
    {_GR_METHOD_LENGTH,     (gr_funcptr) _gr_fmpz_mod_mpoly_q_length},
    {GR_METHOD_WRITE,       (gr_funcptr) _gr_fmpz_mod_mpoly_q_write},
    {GR_METHOD_ZERO,        (gr_funcptr) _gr_fmpz_mod_mpoly_q_zero},
    {GR_METHOD_ONE,         (gr_funcptr) _gr_fmpz_mod_mpoly_q_one},
    {GR_METHOD_IS_ZERO,     (gr_funcptr) _gr_fmpz_mod_mpoly_q_is_zero},
    {GR_METHOD_IS_ONE,      (gr_funcptr) _gr_fmpz_mod_mpoly_q_is_one},
    {GR_METHOD_GENS,        (gr_funcptr) _gr_fmpz_mod_mpoly_q_gens},
    {GR_METHOD_EQUAL,       (gr_funcptr) _gr_fmpz_mod_mpoly_q_equal},
    {GR_METHOD_SET,         (gr_funcptr) _gr_fmpz_mod_mpoly_q_set},
    {GR_METHOD_SET_UI,      (gr_funcptr) _gr_fmpz_mod_mpoly_q_set_ui},
    {GR_METHOD_SET_SI,      (gr_funcptr) _gr_fmpz_mod_mpoly_q_set_si},
    {GR_METHOD_SET_FMPZ,    (gr_funcptr) _gr_fmpz_mod_mpoly_q_set_fmpz},
    {GR_METHOD_SET_FMPQ,    (gr_funcptr) _gr_fmpz_mod_mpoly_q_set_fmpq},
    {GR_METHOD_SET_STR,     (gr_funcptr) gr_generic_set_str_balance_additions},
    {GR_METHOD_NEG,         (gr_funcptr) _gr_fmpz_mod_mpoly_q_neg},
    {GR_METHOD_ADD,         (gr_funcptr) _gr_fmpz_mod_mpoly_q_add},
    {GR_METHOD_ADD_SI,      (gr_funcptr) _gr_fmpz_mod_mpoly_q_add_si},

    {GR_METHOD_ADD_FMPZ,    (gr_funcptr) _gr_fmpz_mod_mpoly_q_add_fmpz},
    {GR_METHOD_ADD_FMPQ,    (gr_funcptr) _gr_fmpz_mod_mpoly_q_add_fmpq},
    {GR_METHOD_SUB,         (gr_funcptr) _gr_fmpz_mod_mpoly_q_sub},
    {GR_METHOD_SUB_SI,      (gr_funcptr) _gr_fmpz_mod_mpoly_q_sub_si},

    {GR_METHOD_SUB_FMPZ,    (gr_funcptr) _gr_fmpz_mod_mpoly_q_sub_fmpz},
    {GR_METHOD_SUB_FMPQ,    (gr_funcptr) _gr_fmpz_mod_mpoly_q_sub_fmpq},
    {GR_METHOD_MUL,         (gr_funcptr) _gr_fmpz_mod_mpoly_q_mul},
    {GR_METHOD_MUL_SI,      (gr_funcptr) _gr_fmpz_mod_mpoly_q_mul_si},

    {GR_METHOD_MUL_FMPZ,    (gr_funcptr) _gr_fmpz_mod_mpoly_q_mul_fmpz},
    {GR_METHOD_MUL_FMPQ,    (gr_funcptr) _gr_fmpz_mod_mpoly_q_mul_fmpq},
    {GR_METHOD_DIV,         (gr_funcptr) _gr_fmpz_mod_mpoly_q_div},
    {GR_METHOD_DIV_SI,      (gr_funcptr) _gr_fmpz_mod_mpoly_q_div_si},
 
    {GR_METHOD_DIV_FMPZ,    (gr_funcptr) _gr_fmpz_mod_mpoly_q_div_fmpz},
    {GR_METHOD_DIV_FMPQ,    (gr_funcptr) _gr_fmpz_mod_mpoly_q_div_fmpq},
    {GR_METHOD_DIVEXACT,         (gr_funcptr) _gr_fmpz_mod_mpoly_q_div},
    {GR_METHOD_DIVEXACT_SI,      (gr_funcptr) _gr_fmpz_mod_mpoly_q_div_si},
 
    {GR_METHOD_DIVEXACT_FMPZ,    (gr_funcptr) _gr_fmpz_mod_mpoly_q_div_fmpz},
    {GR_METHOD_DIVEXACT_FMPQ,    (gr_funcptr) _gr_fmpz_mod_mpoly_q_div_fmpq},
    {GR_METHOD_INV,             (gr_funcptr) _gr_fmpz_mod_mpoly_q_inv},
    {GR_METHOD_IS_INVERTIBLE,   (gr_funcptr) _gr_fmpz_mod_mpoly_q_is_invertible},
    {GR_METHOD_POW_UI,      (gr_funcptr) _gr_fmpz_mod_mpoly_q_pow_ui},
    {GR_METHOD_POW_FMPZ,    (gr_funcptr) _gr_fmpz_mod_mpoly_q_pow_fmpz},
    {GR_METHOD_NUMERATOR,      (gr_funcptr) _gr_fmpz_mod_mpoly_q_numerator},
    {GR_METHOD_DENOMINATOR,      (gr_funcptr) _gr_fmpz_mod_mpoly_q_denominator},
    {0,                     (gr_funcptr) NULL},
};



void
gr_ctx_init_fmpz_mod_mpoly_q(gr_ctx_t ctx, slong nvars, const ordering_t ord, const fmpz_t mod)
{
    ctx->which_ring = GR_CTX_FMPZ_MOD_MPOLY_Q;
    ctx->sizeof_elem = sizeof(fmpz_mod_mpoly_q_struct);
    GR_CTX_DATA_AS_PTR(ctx) = flint_malloc(sizeof(_gr_fmpz_mod_mpoly_ctx_t));
    ctx->size_limit = WORD_MAX;

    fmpz_mod_mpoly_ctx_init(MPOLYNOMIAL_MCTX(ctx), nvars, ord, mod);
    MPOLYNOMIAL_CTX(ctx)->vars = NULL;

    ctx->methods = _gr_fmpz_mod_mpoly_q_methods;

    if (!_gr_fmpz_mod_mpoly_q_methods_initialized)
    {
        gr_method_tab_init(_gr_fmpz_mod_mpoly_q_methods, _gr_fmpz_mod_mpoly_q_methods_input);
        _gr_fmpz_mod_mpoly_q_methods_initialized = 1;
    }
}
