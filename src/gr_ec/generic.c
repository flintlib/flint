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
    The generic gr interface for the group of points of a curve.

    A curve is not a ring: there is no multiplication of points, and no
    integer can be cast to a point. Those operations are mapped to
    gr_not_in_domain (GR_DOMAIN, "no such thing here") rather than left at
    the default gr_not_implemented (GR_UNABLE, "cannot compute it"), so
    that callers can tell the two apart.

    What the interface does provide is the abelian group structure, with
    gr_mul_ui, gr_mul_si, gr_mul_fmpz, gr_mul_two and gr_mul_2exp_si /
    gr_mul_2exp_fmpz acting as the Z-module scalar multiplication. That
    agrees with the documented meaning of those methods, since x*(n*1) is
    n copies of x added in any ring.

    The 2exp methods reject a negative exponent with GR_DOMAIN, since
    halving a point is not single valued.
*/

/* Context predicates */

static truth_t
_gr_ec_ctx_is_finite(gr_ec_ctx_t ctx)
{
    return gr_ctx_is_finite(GR_EC_ELEM_CTX(ctx));
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

/*
    gr_randtest is required to succeed, but finding a point on a curve can
    fail: over characteristic 2 there is no x-coordinate lift, and over an
    infinite ring the search may simply run out of attempts. Fall back to
    the point at infinity, which is always available.
*/

#define GR_EC_RANDTEST_WRAPPER(name, kind) \
static int \
name(kind ## _t res, flint_rand_t state, gr_ec_ctx_t ctx) \
{ \
    int status = kind ## _randtest(res, state, ctx); \
    if (status != GR_SUCCESS) \
        status = kind ## _zero(res, ctx); \
    return status; \
}

GR_EC_RANDTEST_WRAPPER(_gr_ec_point_randtest, gr_ec_point)
GR_EC_RANDTEST_WRAPPER(_gr_ec_aff_point_randtest, gr_ec_aff_point)
GR_EC_RANDTEST_WRAPPER(_gr_ec_jac_point_randtest, gr_ec_jac_point)

/* Entries shared by all three representations. */

#define GR_EC_CTX_METHODS \
    {GR_METHOD_CTX_WRITE,           (gr_funcptr) gr_ec_ctx_write}, \
    {GR_METHOD_CTX_CLEAR,           (gr_funcptr) gr_ec_ctx_clear}, \
    {GR_METHOD_CTX_IS_RING,         (gr_funcptr) gr_generic_ctx_predicate_false}, \
    {GR_METHOD_CTX_IS_COMMUTATIVE_RING, (gr_funcptr) gr_generic_ctx_predicate_false}, \
    {GR_METHOD_CTX_IS_INTEGRAL_DOMAIN,  (gr_funcptr) gr_generic_ctx_predicate_false}, \
    {GR_METHOD_CTX_IS_UNIQUE_FACTORIZATION_DOMAIN, (gr_funcptr) gr_generic_ctx_predicate_false}, \
    {GR_METHOD_CTX_IS_FIELD,        (gr_funcptr) gr_generic_ctx_predicate_false}, \
    {GR_METHOD_CTX_IS_ALGEBRAICALLY_CLOSED, (gr_funcptr) gr_generic_ctx_predicate_false}, \
    {GR_METHOD_CTX_IS_ORDERED_RING, (gr_funcptr) gr_generic_ctx_predicate_false}, \
    {GR_METHOD_CTX_IS_MULTIPLICATIVE_GROUP, (gr_funcptr) gr_generic_ctx_predicate_false}, \
    {GR_METHOD_CTX_IS_RATIONAL_VECTOR_SPACE, (gr_funcptr) gr_generic_ctx_predicate_false}, \
    {GR_METHOD_CTX_IS_REAL_VECTOR_SPACE, (gr_funcptr) gr_generic_ctx_predicate_false}, \
    {GR_METHOD_CTX_IS_COMPLEX_VECTOR_SPACE, (gr_funcptr) gr_generic_ctx_predicate_false}, \
    /* infinity has arbitrary coordinates, and (X:Y:Z) is a class */ \
    {GR_METHOD_CTX_IS_CANONICAL,    (gr_funcptr) gr_generic_ctx_predicate_false}, \
    {GR_METHOD_CTX_IS_FINITE,       (gr_funcptr) _gr_ec_ctx_is_finite}, \
    {GR_METHOD_CTX_IS_EXACT,        (gr_funcptr) _gr_ec_ctx_is_exact}, \
    {GR_METHOD_CTX_IS_THREADSAFE,   (gr_funcptr) _gr_ec_ctx_is_threadsafe}, \
    {GR_METHOD_CTX_BASE,            (gr_funcptr) _gr_ec_ctx_base},

/*
    Undefined on a curve, as opposed to merely unimplemented: there is no
    multiplication, no identity for it, and no way to read an integer or a
    rational as a point. Division by n is a genuine operation on a curve,
    but it needs division polynomials and has several answers, so it does
    not belong in a slot that means "multiply by a rational".
*/

#define GR_EC_NOT_IN_DOMAIN_METHODS \
    {GR_METHOD_ONE,         (gr_funcptr) gr_not_in_domain}, \
    {GR_METHOD_NEG_ONE,     (gr_funcptr) gr_not_in_domain}, \
    {GR_METHOD_SET_UI,      (gr_funcptr) gr_not_in_domain}, \
    {GR_METHOD_SET_SI,      (gr_funcptr) gr_not_in_domain}, \
    {GR_METHOD_SET_FMPZ,    (gr_funcptr) gr_not_in_domain}, \
    {GR_METHOD_SET_FMPQ,    (gr_funcptr) gr_not_in_domain}, \
    {GR_METHOD_SET_STR,     (gr_funcptr) gr_not_in_domain}, \
    {GR_METHOD_MUL,         (gr_funcptr) gr_not_in_domain}, \
    {GR_METHOD_SQR,         (gr_funcptr) gr_not_in_domain}, \
    {GR_METHOD_MUL_FMPQ,    (gr_funcptr) gr_not_in_domain}, \
    {GR_METHOD_DIV,         (gr_funcptr) gr_not_in_domain}, \
    {GR_METHOD_INV,         (gr_funcptr) gr_not_in_domain}, \
    {GR_METHOD_POW_UI,      (gr_funcptr) gr_not_in_domain}, \
    {GR_METHOD_POW_SI,      (gr_funcptr) gr_not_in_domain}, \
    {GR_METHOD_POW_FMPZ,    (gr_funcptr) gr_not_in_domain},

/* The group operations, one block per representation. */

#define GR_EC_GROUP_METHODS(kind, randtest) \
    {GR_METHOD_INIT,        (gr_funcptr) kind ## _init}, \
    {GR_METHOD_CLEAR,       (gr_funcptr) kind ## _clear}, \
    {GR_METHOD_SWAP,        (gr_funcptr) kind ## _swap}, \
    {GR_METHOD_RANDTEST,    (gr_funcptr) randtest}, \
    {GR_METHOD_WRITE,       (gr_funcptr) kind ## _write}, \
    {GR_METHOD_ZERO,        (gr_funcptr) kind ## _zero}, \
    {GR_METHOD_IS_ZERO,     (gr_funcptr) kind ## _is_inf}, \
    {GR_METHOD_EQUAL,       (gr_funcptr) kind ## _equal}, \
    {GR_METHOD_SET,         (gr_funcptr) kind ## _set}, \
    {GR_METHOD_NEG,         (gr_funcptr) kind ## _neg}, \
    {GR_METHOD_ADD,         (gr_funcptr) kind ## _add}, \
    {GR_METHOD_SUB,         (gr_funcptr) kind ## _sub}, \
    {GR_METHOD_MUL_TWO,     (gr_funcptr) kind ## _dbl}, \
    {GR_METHOD_MUL_UI,      (gr_funcptr) kind ## _mul_ui}, \
    {GR_METHOD_MUL_SI,      (gr_funcptr) kind ## _mul_si}, \
    {GR_METHOD_MUL_FMPZ,    (gr_funcptr) kind ## _mul_fmpz}, \
    {GR_METHOD_MUL_2EXP_SI,   (gr_funcptr) kind ## _mul_2exp_si}, \
    {GR_METHOD_MUL_2EXP_FMPZ, (gr_funcptr) kind ## _mul_2exp_fmpz},

static int _gr_ec_methods_initialized[GR_EC_NUM_REPRS] = { 0, 0, 0 };

static gr_static_method_table _gr_ec_methods[GR_EC_NUM_REPRS];

static gr_method_tab_input _gr_ec_projective_methods_input[] =
{
    GR_EC_CTX_METHODS
    GR_EC_GROUP_METHODS(gr_ec_point, _gr_ec_point_randtest)
    GR_EC_NOT_IN_DOMAIN_METHODS
    {0,                     (gr_funcptr) NULL},
};

static gr_method_tab_input _gr_ec_affine_methods_input[] =
{
    GR_EC_CTX_METHODS
    GR_EC_GROUP_METHODS(gr_ec_aff_point, _gr_ec_aff_point_randtest)
    GR_EC_NOT_IN_DOMAIN_METHODS
    {0,                     (gr_funcptr) NULL},
};

static gr_method_tab_input _gr_ec_jacobian_methods_input[] =
{
    GR_EC_CTX_METHODS
    GR_EC_GROUP_METHODS(gr_ec_jac_point, _gr_ec_jac_point_randtest)
    GR_EC_NOT_IN_DOMAIN_METHODS
    {0,                     (gr_funcptr) NULL},
};

static gr_method_tab_input * _gr_ec_methods_input[GR_EC_NUM_REPRS] =
{
    _gr_ec_projective_methods_input,
    _gr_ec_affine_methods_input,
    _gr_ec_jacobian_methods_input,
};

static const slong _gr_ec_sizeof_elem[GR_EC_NUM_REPRS] =
{
    sizeof(gr_ec_point_struct),
    sizeof(gr_ec_aff_point_struct),
    sizeof(gr_ec_jac_point_struct),
};

void
_gr_ec_ctx_init_methods(gr_ec_ctx_t ctx, gr_ec_repr_t repr)
{
    GR_EC_CTX(ctx)->repr = repr;

    ctx->sizeof_elem = _gr_ec_sizeof_elem[repr];
    ctx->methods = _gr_ec_methods[repr];

    if (!_gr_ec_methods_initialized[repr])
    {
        gr_method_tab_init(_gr_ec_methods[repr], _gr_ec_methods_input[repr]);
        _gr_ec_methods_initialized[repr] = 1;
    }
}

int
gr_ec_ctx_set_repr(gr_ec_ctx_t ctx, gr_ec_repr_t repr)
{
    if (repr < 0 || repr >= GR_EC_NUM_REPRS)
        return GR_DOMAIN;

    /* Affine arithmetic inverts coordinates. */
    if (repr == GR_EC_REPR_AFFINE
            && gr_ctx_is_field(GR_EC_ELEM_CTX(ctx)) == T_FALSE)
        return GR_DOMAIN;

    _gr_ec_ctx_init_methods(ctx, repr);

    return GR_SUCCESS;
}

int
gr_ctx_init_gr_ec(gr_ctx_t ctx, gr_ctx_t base_ring, gr_srcptr a1, gr_srcptr a2,
        gr_srcptr a3, gr_srcptr a4, gr_srcptr a6, gr_ec_repr_t repr)
{
    int status;

    if (repr < 0 || repr >= GR_EC_NUM_REPRS)
        return GR_DOMAIN;

    if (repr == GR_EC_REPR_AFFINE && gr_ctx_is_field(base_ring) == T_FALSE)
        return GR_DOMAIN;

    status = gr_ec_ctx_init(ctx, base_ring, a1, a2, a3, a4, a6);

    if (status == GR_SUCCESS)
        _gr_ec_ctx_init_methods(ctx, repr);

    return status;
}
