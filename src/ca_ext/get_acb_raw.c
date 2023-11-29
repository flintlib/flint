/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_hypgeom.h"
#include "ca_ext.h"

#define ARB_CONST(f) \
    f(acb_realref(res), prec); \
    arb_zero(acb_imagref(res)); \
    break;

#define ACB_UNARY(f) \
    ca_get_acb_raw(res, CA_EXT_FUNC_ARGS(x), prec, ctx); \
    f(res, res, prec); \
    break;

#define ACB_UNARY_NOPREC(f) \
    ca_get_acb_raw(res, CA_EXT_FUNC_ARGS(x), prec, ctx); \
    f(res, res); \
    break;

#define ACB_UNARY_REAL(f) \
    ca_get_acb_raw(res, CA_EXT_FUNC_ARGS(x), prec, ctx); \
    f(acb_realref(res), res, prec); \
    arb_zero(acb_imagref(res)); \
    break;

#define ACB_UNARY_REAL_NOPREC(f) \
    ca_get_acb_raw(res, CA_EXT_FUNC_ARGS(x), prec, ctx); \
    f(acb_realref(res), res); \
    arb_zero(acb_imagref(res)); \
    break;

#define ACB_UNARY_REAL_REAL(f) \
    ca_get_acb_raw(res, CA_EXT_FUNC_ARGS(x), prec, ctx); \
    f(acb_realref(res), acb_realref(res), prec); \
    arb_zero(acb_imagref(res)); \
    break;

#define ACB_BINARY(f) \
    { \
        acb_t _t; \
        acb_init(_t); \
        ca_get_acb_raw(res, CA_EXT_FUNC_ARGS(x), prec, ctx); \
        ca_get_acb_raw(_t, CA_EXT_FUNC_ARGS(x) + 1, prec, ctx); \
        f(res, res, _t, prec); \
        acb_clear(_t); \
    } \
    break;

void
ca_ext_get_acb_raw(acb_t res, ca_ext_t x, slong prec, ca_ctx_t ctx)
{
    if (CA_EXT_HEAD(x) == CA_QQBar)
    {
        qqbar_cache_enclosure(CA_EXT_QQBAR(x), prec);
        qqbar_get_acb(res, CA_EXT_QQBAR(x), prec);
        return;
    }

    if (prec <= CA_EXT_FUNC_PREC(x))
    {
        acb_set(res, CA_EXT_FUNC_ENCLOSURE(x));
        return;
    }

    switch (CA_EXT_HEAD(x))
    {
        /* Arithmetic */
        case CA_Neg: ACB_UNARY_NOPREC(acb_neg)
        case CA_Add: ACB_BINARY(acb_add)
        case CA_Sub: ACB_BINARY(acb_sub)
        case CA_Mul: ACB_BINARY(acb_mul)
        case CA_Div: ACB_BINARY(acb_div)
        /* Roots */
        case CA_Sqrt: ACB_UNARY(acb_sqrt)
        /* CA_Cbrt,  not implemented */
        /* CA_Root,  not implemented */
        /* Complex parts */
        case CA_Floor: ACB_UNARY_REAL_REAL(arb_floor)
        case CA_Ceil: ACB_UNARY_REAL_REAL(arb_ceil)
        case CA_Abs:  ACB_UNARY_REAL(acb_abs)
        case CA_Sign: ACB_UNARY(acb_sgn)
        case CA_Re:   ACB_UNARY_REAL_NOPREC(acb_get_real)
        case CA_Im:   ACB_UNARY_REAL_NOPREC(acb_get_imag)
        case CA_Arg:  ACB_UNARY_REAL(acb_arg)
        case CA_Conjugate:  ACB_UNARY_NOPREC(acb_conj)
        /* Elementary constants */
        case CA_Pi:   ARB_CONST(arb_const_pi)
        /* Elementary functions */
        case CA_Exp:   ACB_UNARY(acb_exp)
        case CA_Log:   ACB_UNARY(acb_log)
        case CA_Pow:   ACB_BINARY(acb_pow)
        case CA_Cos:   ACB_UNARY(acb_cos)
        case CA_Sin:   ACB_UNARY(acb_sin)
        case CA_Tan:   ACB_UNARY(acb_tan)
        case CA_Cosh:  ACB_UNARY(acb_cosh)
        case CA_Sinh:  ACB_UNARY(acb_sinh)
        case CA_Tanh:  ACB_UNARY(acb_tanh)
        case CA_Atan:   ACB_UNARY(acb_atan)
        case CA_Acos:   ACB_UNARY(acb_acos)
        case CA_Asin:   ACB_UNARY(acb_asin)
        case CA_Atanh:  ACB_UNARY(acb_atanh)
        case CA_Acosh:  ACB_UNARY(acb_acosh)
        case CA_Asinh:  ACB_UNARY(acb_asinh)
        /* Euler's constant */
        case CA_Euler: ARB_CONST(arb_const_euler)
        /* Gamma and related functions */
        case CA_Gamma:        ACB_UNARY(acb_gamma)
        case CA_LogGamma:     ACB_UNARY(acb_lgamma)
        case CA_Psi:          ACB_UNARY(acb_digamma)
        case CA_Erf:          ACB_UNARY(acb_hypgeom_erf)
        case CA_Erfc:         ACB_UNARY(acb_hypgeom_erfc)
        case CA_Erfi:         ACB_UNARY(acb_hypgeom_erfi)
        case CA_RiemannZeta:  ACB_UNARY(acb_zeta)
        case CA_HurwitzZeta:  ACB_BINARY(acb_hurwitz_zeta)
        default:
            flint_throw(FLINT_ERROR, "ca_ext_get_acb_raw: unknown function\n");
    }

    acb_set(CA_EXT_FUNC_ENCLOSURE(x), res);
    CA_EXT_FUNC_PREC(x) = prec;
}

