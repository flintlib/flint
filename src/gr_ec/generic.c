/*
    Copyright (C) 2026 Maël Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr.h"
#include "gr_generic.h"
#include "gr_ec.h"
#include "impl.h"

/*
    The generic interface to the group of points, whose elements are
    Jacobian points. The ring operations do not exist on a curve and return
    GR_DOMAIN; the integer scalar multiplications are the Z-module action.
*/

/* E(R) is finite over a finite R, but may also be finite over an infinite
   one (y^2 = x^3 - x has four points over Q) */
static truth_t
_gr_ec_ctx_is_finite(gr_ec_ctx_t ctx)
{
    return (gr_ctx_is_finite(GR_EC_ELEM_CTX(ctx)) == T_TRUE) ? T_TRUE : T_UNKNOWN;
}

static truth_t
_gr_ec_ctx_is_exact(gr_ec_ctx_t ctx)
{
    return gr_ctx_is_exact(GR_EC_ELEM_CTX(ctx));
}

static truth_t
_gr_ec_ctx_is_threadsafe(gr_ec_ctx_t ctx)
{
    return gr_ctx_is_threadsafe(GR_EC_ELEM_CTX(ctx));
}

static gr_ctx_ptr
_gr_ec_ctx_base(gr_ec_ctx_t ctx)
{
    return GR_EC_ELEM_CTX(ctx);
}

/* reads the cached count but never writes it: no gr method modifies ctx */
static int
_gr_ec_ctx_cardinality_fmpz(fmpz_t res, gr_ec_ctx_t ctx)
{
    if (!fmpz_is_zero(&GR_EC_CTX(ctx)->cardinality))
    {
        fmpz_set(res, &GR_EC_CTX(ctx)->cardinality);
        return GR_SUCCESS;
    }

    return _gr_ec_ctx_cardinality(res, ctx);
}

/* gr_randtest must succeed, but finding a point can fail (for instance in
   characteristic 2 over an infinite ring), so fall back to the identity */
static int
_gr_ec_jac_point_randtest(gr_ec_jac_point_t res, flint_rand_t state,
        gr_ec_ctx_t ctx)
{
    if (gr_ec_jac_point_randtest(res, state, ctx) != GR_SUCCESS)
        return gr_ec_jac_point_zero(res, ctx);

    return GR_SUCCESS;
}

static int _gr_ec_methods_initialized = 0;

static gr_static_method_table _gr_ec_methods;

static gr_method_tab_input _gr_ec_methods_input[] =
{
    {GR_METHOD_CTX_WRITE,           (gr_funcptr) _gr_ec_ctx_write},
    {GR_METHOD_CTX_CLEAR,           (gr_funcptr) gr_ec_ctx_clear},
    {GR_METHOD_CTX_IS_RING,         (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_COMMUTATIVE_RING, (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_INTEGRAL_DOMAIN,  (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_UNIQUE_FACTORIZATION_DOMAIN, (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_FIELD,        (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_ALGEBRAICALLY_CLOSED, (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_ORDERED_RING, (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_MULTIPLICATIVE_GROUP, (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_RATIONAL_VECTOR_SPACE, (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_REAL_VECTOR_SPACE, (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_COMPLEX_VECTOR_SPACE, (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_CANONICAL,    (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_FINITE,       (gr_funcptr) _gr_ec_ctx_is_finite},
    {GR_METHOD_CTX_IS_EXACT,        (gr_funcptr) _gr_ec_ctx_is_exact},
    {GR_METHOD_CTX_IS_THREADSAFE,   (gr_funcptr) _gr_ec_ctx_is_threadsafe},
    {GR_METHOD_CTX_BASE,            (gr_funcptr) _gr_ec_ctx_base},
    {GR_METHOD_CTX_CARDINALITY_FMPZ, (gr_funcptr) _gr_ec_ctx_cardinality_fmpz},

    {GR_METHOD_INIT,                (gr_funcptr) gr_ec_jac_point_init},
    {GR_METHOD_CLEAR,               (gr_funcptr) gr_ec_jac_point_clear},
    {GR_METHOD_SWAP,                (gr_funcptr) gr_ec_jac_point_swap},
    {GR_METHOD_RANDTEST,            (gr_funcptr) _gr_ec_jac_point_randtest},
    {GR_METHOD_WRITE,               (gr_funcptr) gr_ec_jac_point_write},
    {GR_METHOD_ZERO,                (gr_funcptr) gr_ec_jac_point_zero},
    {GR_METHOD_IS_ZERO,             (gr_funcptr) gr_ec_jac_point_is_inf},
    {GR_METHOD_EQUAL,               (gr_funcptr) gr_ec_jac_point_equal},
    {GR_METHOD_SET,                 (gr_funcptr) gr_ec_jac_point_set},
    {GR_METHOD_NEG,                 (gr_funcptr) gr_ec_jac_point_neg},
    {GR_METHOD_ADD,                 (gr_funcptr) gr_ec_jac_point_add},
    {GR_METHOD_SUB,                 (gr_funcptr) gr_ec_jac_point_sub},
    {GR_METHOD_MUL_TWO,             (gr_funcptr) gr_ec_jac_point_dbl},
    {GR_METHOD_MUL_UI,              (gr_funcptr) gr_ec_jac_point_mul_ui},
    {GR_METHOD_MUL_SI,              (gr_funcptr) gr_ec_jac_point_mul_si},
    {GR_METHOD_MUL_FMPZ,            (gr_funcptr) gr_ec_jac_point_mul_fmpz},
    {GR_METHOD_MUL_2EXP_SI,         (gr_funcptr) gr_ec_jac_point_mul_2exp_si},
    {GR_METHOD_MUL_2EXP_FMPZ,       (gr_funcptr) gr_ec_jac_point_mul_2exp_fmpz},

    {GR_METHOD_ONE,                 (gr_funcptr) gr_not_in_domain},
    {GR_METHOD_NEG_ONE,             (gr_funcptr) gr_not_in_domain},
    {GR_METHOD_SET_UI,              (gr_funcptr) gr_not_in_domain},
    {GR_METHOD_SET_SI,              (gr_funcptr) gr_not_in_domain},
    {GR_METHOD_SET_FMPZ,            (gr_funcptr) gr_not_in_domain},
    {GR_METHOD_SET_FMPQ,            (gr_funcptr) gr_not_in_domain},
    {GR_METHOD_SET_STR,             (gr_funcptr) gr_not_in_domain},
    {GR_METHOD_MUL,                 (gr_funcptr) gr_not_in_domain},
    {GR_METHOD_SQR,                 (gr_funcptr) gr_not_in_domain},
    {GR_METHOD_DIV,                 (gr_funcptr) gr_not_in_domain},
    {GR_METHOD_INV,                 (gr_funcptr) gr_not_in_domain},
    {GR_METHOD_POW_UI,              (gr_funcptr) gr_not_in_domain},
    {GR_METHOD_POW_SI,              (gr_funcptr) gr_not_in_domain},
    {GR_METHOD_POW_FMPZ,            (gr_funcptr) gr_not_in_domain},
    {0,                             (gr_funcptr) NULL},
};

void
_gr_ec_ctx_init_methods(gr_ec_ctx_t ctx)
{
    ctx->which_ring = GR_CTX_GR_EC;
    ctx->sizeof_elem = sizeof(gr_ec_jac_point_struct);
    ctx->size_limit = WORD_MAX;
    ctx->methods = _gr_ec_methods;

    if (!_gr_ec_methods_initialized)
    {
        gr_method_tab_init(_gr_ec_methods, _gr_ec_methods_input);
        _gr_ec_methods_initialized = 1;
    }
}
