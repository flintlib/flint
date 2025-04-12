/*
    Copyright (C) 2025 Ricardo Buring

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Ore polynomials over generic rings */

#include <stdlib.h>
#include <string.h>
#include "gr_vec.h"
#include "gr_ore_poly.h"
#include "gr_generic.h"

static const char * default_var = "D";

int gr_ore_poly_ctx_write(gr_stream_t out, gr_ore_poly_ctx_t ctx)
{
    gr_stream_write(out, "Ring of Ore polynomials over ");
    gr_ctx_write(out, GR_ORE_POLY_ELEM_CTX(ctx));
    return GR_SUCCESS;
}

int _gr_ore_poly_ctx_set_gen_name(gr_ctx_t ctx, const char * s)
{
    slong len;
    len = strlen(s);

    if (GR_ORE_POLY_CTX(ctx)->var == default_var)
        GR_ORE_POLY_CTX(ctx)->var = NULL;

    GR_ORE_POLY_CTX(ctx)->var = flint_realloc(GR_ORE_POLY_CTX(ctx)->var, len + 1);
    memcpy(GR_ORE_POLY_CTX(ctx)->var, s, len + 1);
    return GR_SUCCESS;
}

int _gr_ore_poly_ctx_set_gen_names(gr_ctx_t ctx, const char ** s)
{
    return _gr_ore_poly_ctx_set_gen_name(ctx, s[0]);
}

void
gr_ore_poly_ctx_clear(gr_ore_poly_ctx_t ctx)
{
    if (GR_ORE_POLY_CTX(ctx)->var != default_var)
    {
        flint_free(GR_ORE_POLY_CTX(ctx)->var);
    }
}

truth_t
gr_ore_poly_ctx_is_ring(gr_ore_poly_ctx_t ctx)
{
    return gr_ctx_is_ring(GR_ORE_POLY_ELEM_CTX(ctx));
}

truth_t
gr_ore_poly_ctx_is_zero_ring(gr_ore_poly_ctx_t ctx)
{
    return gr_ctx_is_zero_ring(GR_ORE_POLY_ELEM_CTX(ctx));
}

truth_t
gr_ore_poly_ctx_is_commutative_ring(gr_ore_poly_ctx_t ctx)
{
    return T_UNKNOWN;
}

truth_t
gr_ore_poly_ctx_is_integral_domain(gr_ore_poly_ctx_t ctx)
{
    return T_UNKNOWN;
}

truth_t
gr_ore_poly_ctx_is_unique_factorization_domain(gr_ore_poly_ctx_t ctx)
{
    return T_UNKNOWN;
}

truth_t
gr_ore_poly_ctx_is_threadsafe(gr_ore_poly_ctx_t ctx)
{
    return gr_ctx_is_threadsafe(GR_ORE_POLY_ELEM_CTX(ctx));
}

int
ore_poly_write(gr_stream_t out, gr_ore_poly_t poly, gr_ctx_t ctx)
{
    /* todo */
    if (poly->length == 0)
    {
        gr_stream_write(out, "0");
        return GR_SUCCESS;
    }

    return gr_ore_poly_write(out, poly, ctx);
}

int
gr_ore_poly_i(gr_ore_poly_t res, gr_ore_poly_ctx_t ctx)
{
    int status;
    gr_ore_poly_fit_length(res, 1, ctx);
    _gr_ore_poly_set_length(res, 1, ctx);
    status = gr_i(res->coeffs, GR_ORE_POLY_ELEM_CTX(ctx));
    _gr_ore_poly_normalise(res, ctx);
    return status;
}

int
gr_ore_poly_gens_recursive(gr_vec_t vec, gr_ore_poly_ctx_t ctx)
{
    int status;
    gr_vec_t vec1;
    slong i, n;

    /* Get generators of the element ring */
    gr_vec_init(vec1, 0, GR_ORE_POLY_ELEM_CTX(ctx));
    status = gr_gens_recursive(vec1, GR_ORE_POLY_ELEM_CTX(ctx));
    n = vec1->length;

    gr_vec_set_length(vec, n + 1, ctx);

    /* Promote to Ore polynomials */
    for (i = 0; i < n; i++)
        status |= gr_ore_poly_set_other(gr_vec_entry_ptr(vec, i, ctx),
                gr_vec_entry_srcptr(vec1, i, GR_ORE_POLY_ELEM_CTX(ctx)),
                GR_ORE_POLY_ELEM_CTX(ctx),
                ctx);

    status |= gr_ore_poly_gen(gr_vec_entry_ptr(vec, n, ctx), ctx);

    gr_vec_clear(vec1, GR_ORE_POLY_ELEM_CTX(ctx));

    return status;
}

int _gr_ore_poly_methods_initialized = 0;

gr_static_method_table _gr_ore_poly_methods;

gr_method_tab_input _gr_ore_poly_methods_input[] =
{
    {GR_METHOD_CTX_WRITE,   (gr_funcptr) gr_ore_poly_ctx_write},
    {GR_METHOD_CTX_CLEAR,   (gr_funcptr) gr_ore_poly_ctx_clear},

    {GR_METHOD_CTX_IS_RING,     (gr_funcptr) gr_ore_poly_ctx_is_ring},
    {GR_METHOD_CTX_IS_ZERO_RING,     (gr_funcptr) gr_ore_poly_ctx_is_zero_ring},
    {GR_METHOD_CTX_IS_COMMUTATIVE_RING, (gr_funcptr) gr_ore_poly_ctx_is_commutative_ring},
    {GR_METHOD_CTX_IS_INTEGRAL_DOMAIN,  (gr_funcptr) gr_ore_poly_ctx_is_integral_domain},
    {GR_METHOD_CTX_IS_UNIQUE_FACTORIZATION_DOMAIN,  (gr_funcptr) gr_ore_poly_ctx_is_unique_factorization_domain},
    {GR_METHOD_CTX_IS_FIELD,            (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_THREADSAFE,       (gr_funcptr) gr_ore_poly_ctx_is_threadsafe},
    {GR_METHOD_CTX_SET_GEN_NAME,        (gr_funcptr) _gr_ore_poly_ctx_set_gen_name},
    {GR_METHOD_CTX_SET_GEN_NAMES,       (gr_funcptr) _gr_ore_poly_ctx_set_gen_names},

    {GR_METHOD_INIT,        (gr_funcptr) gr_ore_poly_init},
    {GR_METHOD_CLEAR,       (gr_funcptr) gr_ore_poly_clear},
    {GR_METHOD_SWAP,        (gr_funcptr) gr_ore_poly_swap},
    {GR_METHOD_SET_SHALLOW, (gr_funcptr) gr_ore_poly_set_shallow},
    {GR_METHOD_RANDTEST,    (gr_funcptr) _gr_ore_poly_randtest_default},
    {GR_METHOD_WRITE,       (gr_funcptr) ore_poly_write},
    {GR_METHOD_ZERO,        (gr_funcptr) gr_ore_poly_zero},
    {GR_METHOD_ONE,         (gr_funcptr) gr_ore_poly_one},
    {GR_METHOD_NEG_ONE,     (gr_funcptr) gr_ore_poly_neg_one},
    {GR_METHOD_I,           (gr_funcptr) gr_ore_poly_i},

    {GR_METHOD_GEN,            (gr_funcptr) gr_ore_poly_gen},
    {GR_METHOD_GENS,           (gr_funcptr) gr_generic_gens_single},
    {GR_METHOD_GENS_RECURSIVE, (gr_funcptr) gr_ore_poly_gens_recursive},

/*
    {GR_METHOD_IS_ZERO,     (gr_funcptr) gr_ore_poly_is_zero},
    {GR_METHOD_IS_ONE,      (gr_funcptr) gr_ore_poly_is_one},
    {GR_METHOD_IS_NEG_ONE,  (gr_funcptr) gr_ore_poly_is_neg_one},
*/
    {GR_METHOD_EQUAL,       (gr_funcptr) gr_ore_poly_equal},
    {GR_METHOD_SET,         (gr_funcptr) gr_ore_poly_set},
    {GR_METHOD_SET_UI,      (gr_funcptr) gr_ore_poly_set_ui},
    {GR_METHOD_SET_SI,      (gr_funcptr) gr_ore_poly_set_si},
    {GR_METHOD_SET_FMPZ,    (gr_funcptr) gr_ore_poly_set_fmpz},
    {GR_METHOD_SET_FMPQ,    (gr_funcptr) gr_ore_poly_set_fmpq},
    {GR_METHOD_SET_OTHER,   (gr_funcptr) gr_ore_poly_set_other},
/*
    {GR_METHOD_SET_INTERVAL_MID_RAD,    (gr_funcptr) gr_ore_poly_set_interval_mid_rad},
*/
    {GR_METHOD_SET_STR,     (gr_funcptr) gr_ore_poly_set_str},
    {GR_METHOD_NEG,         (gr_funcptr) gr_ore_poly_neg},
    {GR_METHOD_ADD_UI,      (gr_funcptr) gr_ore_poly_add_ui},
    {GR_METHOD_ADD_SI,      (gr_funcptr) gr_ore_poly_add_si},
    {GR_METHOD_ADD_FMPZ,    (gr_funcptr) gr_ore_poly_add_fmpz},
    {GR_METHOD_ADD_FMPQ,    (gr_funcptr) gr_ore_poly_add_fmpq},
    {GR_METHOD_ADD,         (gr_funcptr) gr_ore_poly_add},
    {GR_METHOD_ADD_OTHER,   (gr_funcptr) gr_ore_poly_add_other},
    {GR_METHOD_SUB_UI,      (gr_funcptr) gr_ore_poly_sub_ui},
    {GR_METHOD_SUB_SI,      (gr_funcptr) gr_ore_poly_sub_si},
    {GR_METHOD_SUB_FMPZ,    (gr_funcptr) gr_ore_poly_sub_fmpz},
    {GR_METHOD_SUB_FMPQ,    (gr_funcptr) gr_ore_poly_sub_fmpq},
    {GR_METHOD_SUB,         (gr_funcptr) gr_ore_poly_sub},
    {GR_METHOD_SUB_OTHER,   (gr_funcptr) gr_ore_poly_sub_other},
/*
    {GR_METHOD_MUL,         (gr_funcptr) gr_ore_poly_mul},
    {GR_METHOD_MUL_OTHER,   (gr_funcptr) gr_ore_poly_mul_other},
*/
    {GR_METHOD_OTHER_MUL,   (gr_funcptr) gr_ore_poly_other_mul},
    {GR_METHOD_MUL_UI,      (gr_funcptr) gr_ore_poly_mul_ui},
    {GR_METHOD_MUL_SI,      (gr_funcptr) gr_ore_poly_mul_si},
    {GR_METHOD_MUL_FMPZ,    (gr_funcptr) gr_ore_poly_mul_fmpz},
    {GR_METHOD_MUL_FMPQ,    (gr_funcptr) gr_ore_poly_mul_fmpq},
/*
    {GR_METHOD_POW_UI,      (gr_funcptr) gr_ore_poly_pow_ui},
    {GR_METHOD_POW_SI,      (gr_funcptr) gr_ore_poly_pow_si},
    {GR_METHOD_POW_FMPZ,    (gr_funcptr) gr_ore_poly_pow_fmpz},
    {GR_METHOD_DIV,         (gr_funcptr) gr_ore_poly_div},
    {GR_METHOD_INV,         (gr_funcptr) gr_ore_poly_inv},

    {GR_METHOD_EUCLIDEAN_DIV,         (gr_funcptr) gr_ore_poly_euclidean_div},
    {GR_METHOD_EUCLIDEAN_REM,         (gr_funcptr) gr_ore_poly_euclidean_rem},
    {GR_METHOD_EUCLIDEAN_DIVREM,      (gr_funcptr) gr_ore_poly_euclidean_divrem},

    {GR_METHOD_GCD,         (gr_funcptr) gr_ore_poly_gcd},

    {GR_METHOD_FACTOR,      (gr_funcptr) gr_ore_poly_factor},
*/

    {0,                     (gr_funcptr) NULL},
};

void
gr_ore_poly_ctx_init(gr_ore_poly_ctx_t ctx, gr_ctx_t base_ring, slong base_var, const ore_algebra_t which_algebra)
{
    ctx->which_ring = GR_CTX_GR_ORE_POLY;
    ctx->sizeof_elem = sizeof(gr_ore_poly_struct);
    ctx->size_limit = WORD_MAX;

    GR_ORE_POLY_CTX(ctx)->base_ring = (gr_ctx_struct *) base_ring;
    GR_ORE_POLY_CTX(ctx)->degree_limit = WORD_MAX;
    GR_ORE_POLY_CTX(ctx)->var = (char *) default_var;
    GR_ORE_POLY_CTX(ctx)->which_algebra = which_algebra;
    GR_ORE_POLY_CTX(ctx)->base_var = base_var;

    ctx->methods = _gr_ore_poly_methods;

    if (!_gr_ore_poly_methods_initialized)
    {
        gr_method_tab_init(_gr_ore_poly_methods, _gr_ore_poly_methods_input);
        _gr_ore_poly_methods_initialized = 1;
    }
}

void
gr_ore_poly_ctx_init_rand(gr_ore_poly_ctx_t ctx, flint_rand_t state, gr_ctx_t base_ring)
{
    gr_ore_poly_ctx_init(ctx, base_ring, 0, ore_algebra_randtest(state));
}
