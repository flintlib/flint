/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_hypgeom.h"
#include "acb_modular.h"
#include "acb_dirichlet.h"
#include "fexpr.h"
#include "fexpr_builtin.h"

void
_acb_root(acb_t res, const acb_t x, const acb_t y, slong prec)
{
    if (acb_is_int(y) && arf_sgn(arb_midref(acb_realref(y))) > 0 && arf_cmpabs_ui(arb_midref(acb_realref(y)), 1000) <= 0)
    {
        acb_root_ui(res, x, arf_get_si(arb_midref(acb_realref(y)), ARF_RND_DOWN), prec);
    }
    else
    {
        acb_t e;
        acb_init(e);
        acb_inv(e, y, prec);
        acb_pow(res, x, e, prec);
        acb_clear(e);
    }
}

#define REQUIRE_NARGS(n)  if (nargs != n) { success = 0; break; }


#define ACB_FUNCTION_1(acb_func) \
    REQUIRE_NARGS(1) \
    success = fexpr_get_acb_raw(res, arg, prec); \
    acb_func(res, res, prec); \
    break;

#define ACB_FUNCTION_1_ANALYTIC_FLAG(acb_func) \
    REQUIRE_NARGS(1) \
    success = fexpr_get_acb_raw(res, arg, prec); \
    acb_func(res, res, 0, prec); \
    break;

#define ACB_FUNCTION_2(acb_func) \
    REQUIRE_NARGS(2) \
    success = fexpr_get_acb_raw(res, arg, prec); \
    if (success) \
    { \
        acb_init(t); \
        fexpr_view_next(arg); \
        success = fexpr_get_acb_raw(t, arg, prec); \
        if (success) \
            acb_func(res, res, t, prec); \
        acb_clear(t); \
    } \
    break;


int
fexpr_get_acb_raw(acb_t res, const fexpr_t expr, slong prec)
{
    if (fexpr_is_integer(expr))
    {
        fmpz_t c;
        fmpz_init(c);
        fexpr_get_fmpz(c, expr);
        acb_set_round_fmpz(res, c, prec);
        fmpz_clear(c);
        return 1;
    }
    else if (fexpr_is_atom(expr))
    {
        ulong op;

        if (!fexpr_is_any_builtin_symbol(expr))
        {
            acb_indeterminate(res);
            return 0;
        }

        /* todo: clean up these cases */
        op = FEXPR_BUILTIN_ID(expr->data[0]);

        if (op == FEXPR_Pi)
        {
            acb_const_pi(res, prec);
            return 1;
        }

        if (op == FEXPR_NumberI)
        {
            acb_onei(res);
            return 1;
        }

        if (op == FEXPR_NumberE)
        {
            arb_const_e(acb_realref(res), prec);
            arb_zero(acb_imagref(res));
            return 1;
        }

        if (op == FEXPR_Euler)
        {
            arb_const_euler(acb_realref(res), prec);
            arb_zero(acb_imagref(res));
            return 1;
        }

        if (op == FEXPR_CatalanConstant)
        {
            arb_const_catalan(acb_realref(res), prec);
            arb_zero(acb_imagref(res));
            return 1;
        }

        if (op == FEXPR_GoldenRatio)
        {
            arb_sqrt_ui(acb_realref(res), 5, prec);
            arb_add_ui(acb_realref(res), acb_realref(res), 1, prec);
            arb_mul_2exp_si(acb_realref(res), acb_realref(res), -1);
            arb_zero(acb_imagref(res));
            return 1;
        }

        acb_indeterminate(res);
        return 0;
    }
    else
    {
        slong nargs;
        fmpz_t n, m;
        acb_t t, u, v, w;
        fexpr_t func;
        fexpr_t arg;
        fexpr_builtin_symbol op;
        slong i;
        int success = 0;

        nargs = fexpr_nargs(expr);

        fexpr_view_func(func, expr);

        if (!fexpr_is_any_builtin_symbol(func))
        {
            acb_indeterminate(res);
            return 0;
        }

        if (nargs > 0)
            fexpr_view_arg(arg, expr, 0);

        op = FEXPR_BUILTIN_ID(func->data[0]);

        switch (op)
        {
        case FEXPR_Abs:
            REQUIRE_NARGS(1)
            success = fexpr_get_acb_raw(res, arg, prec);
            acb_abs(acb_realref(res), res, prec);
            arb_zero(acb_imagref(res));
            break;

        case FEXPR_Acos:
            ACB_FUNCTION_1(acb_acos)

        case FEXPR_Acosh:
            ACB_FUNCTION_1(acb_acosh)

        case FEXPR_Add:
            if (nargs == 0)
            {
                acb_zero(res);
            }
            else if (nargs == 1)
            {
                success = fexpr_get_acb_raw(res, arg, prec);
            }
            else
            {
                acb_init(t);
                success = fexpr_get_acb_raw(res, arg, prec);
                for (i = 1; success && i < nargs; i++)
                {
                    fexpr_view_next(arg);
                    success = fexpr_get_acb_raw(t, arg, prec);
                    acb_add(res, res, t, prec);
                }
                acb_clear(t);
            }
            break;

        case FEXPR_AiryAi:
            REQUIRE_NARGS(1)
            success = fexpr_get_acb_raw(res, arg, prec);
            acb_hypgeom_airy(res, NULL, NULL, NULL, res, prec);
            break;

        case FEXPR_AiryBi:
            REQUIRE_NARGS(1)
            success = fexpr_get_acb_raw(res, arg, prec);
            acb_hypgeom_airy(NULL, res, NULL, NULL, res, prec);
            break;

        case FEXPR_Arg:
            REQUIRE_NARGS(1)
            success = fexpr_get_acb_raw(res, arg, prec);
            acb_arg(acb_realref(res), res, prec);
            arb_zero(acb_imagref(res));
            break;

        case FEXPR_Asin:
            ACB_FUNCTION_1(acb_asin)

        case FEXPR_Asinh:
            ACB_FUNCTION_1(acb_asinh)

        case FEXPR_Atan:
            ACB_FUNCTION_1(acb_atan)

        case FEXPR_Atanh:
            ACB_FUNCTION_1(acb_atanh)

        case FEXPR_BesselI:
            ACB_FUNCTION_2(acb_hypgeom_bessel_i)

        case FEXPR_BesselJ:
            ACB_FUNCTION_2(acb_hypgeom_bessel_j)

        case FEXPR_BesselK:
            ACB_FUNCTION_2(acb_hypgeom_bessel_k)

        case FEXPR_BesselY:
            ACB_FUNCTION_2(acb_hypgeom_bessel_y)

        case FEXPR_Ceil:
            ACB_FUNCTION_1_ANALYTIC_FLAG(acb_real_ceil)

        case FEXPR_Conjugate:
            REQUIRE_NARGS(1)
            success = fexpr_get_acb_raw(res, arg, prec);
            acb_conj(res, res);
            break;

        case FEXPR_Cos:
            ACB_FUNCTION_1(acb_cos)

        case FEXPR_Cosh:
            ACB_FUNCTION_1(acb_cosh)

        case FEXPR_DedekindEta:
            ACB_FUNCTION_1(acb_modular_eta)

        case FEXPR_Div:
            ACB_FUNCTION_2(acb_div)

        case FEXPR_Erf:
            ACB_FUNCTION_1(acb_hypgeom_erf)

        case FEXPR_Erfc:
            ACB_FUNCTION_1(acb_hypgeom_erfc)

        case FEXPR_Erfi:
            ACB_FUNCTION_1(acb_hypgeom_erfi)

        case FEXPR_Exp:
            ACB_FUNCTION_1(acb_exp)

        case FEXPR_Floor:
            ACB_FUNCTION_1_ANALYTIC_FLAG(acb_real_floor)

        case FEXPR_Gamma:
            ACB_FUNCTION_1(acb_gamma)

        case FEXPR_HurwitzZeta:
            ACB_FUNCTION_2(acb_hurwitz_zeta)

        case FEXPR_Im:
            REQUIRE_NARGS(1)
            success = fexpr_get_acb_raw(res, arg, prec);
            arb_swap(acb_realref(res), acb_imagref(res));
            arb_zero(acb_imagref(res));
            break;

        case FEXPR_JacobiTheta:
            REQUIRE_NARGS(3)
            {
                fmpz_init(n);
                success = fexpr_get_fmpz(n, arg);
                if (success)
                    success = (*n == 1 || *n == 2 || *n == 3 || *n == 4);
                if (success)
                {
                    fexpr_view_next(arg);
                    success = fexpr_get_acb_raw(res, arg, prec);
                    if (success)
                    {
                        acb_init(t);
                        fexpr_view_next(arg);
                        success = fexpr_get_acb_raw(t, arg, prec);
                        if (success)
                        {
                            acb_init(u);
                            acb_init(v);
                            acb_init(w);
                            if (*n == 1)
                                acb_modular_theta(res, t, u, v, res, t, prec);
                            else if (*n == 2)
                                acb_modular_theta(t, res, u, v, res, t, prec);
                            else if (*n == 3)
                                acb_modular_theta(t, u, res, v, res, t, prec);
                            else
                                acb_modular_theta(t, u, v, res, res, t, prec);
                            acb_clear(u);
                            acb_clear(v);
                            acb_clear(w);
                        }
                        acb_clear(t);
                    }
                }
                fmpz_clear(n);
            }
            break;

        case FEXPR_LambertW:
            if (nargs == 1)
            {
                fmpz_init(n);
                success = fexpr_get_acb_raw(res, arg, prec);
                if (success)
                    acb_lambertw(res, res, n, 0, prec);
                fmpz_clear(n);
            }
            else if (nargs == 2)
            {
                fmpz_init(n);
                success = fexpr_get_fmpz(n, arg);
                if (success)
                {
                    fexpr_view_next(arg);
                    success = fexpr_get_acb_raw(res, arg, prec);
                    if (success)
                        acb_lambertw(res, res, n, 0, prec);
                }
                fmpz_clear(n);
            }
            break;

        case FEXPR_Log:
            ACB_FUNCTION_1(acb_log)

        case FEXPR_LogGamma:
            ACB_FUNCTION_1(acb_lgamma)

        case FEXPR_Mul:
            if (nargs == 0)
            {
                acb_one(res);
            }
            else if (nargs == 1)
            {
                success = fexpr_get_acb_raw(res, arg, prec);
            }
            else
            {
                acb_init(t);
                success = fexpr_get_acb_raw(res, arg, prec);
                for (i = 1; success && i < nargs; i++)
                {
                    fexpr_view_next(arg);
                    success = fexpr_get_acb_raw(t, arg, prec);
                    acb_mul(res, res, t, prec);
                }
                acb_clear(t);
            }
            break;

        case FEXPR_Neg:
            REQUIRE_NARGS(1)
            success = fexpr_get_acb_raw(res, arg, prec);
            acb_neg(res, res);
            break;

        case FEXPR_Pos:
            REQUIRE_NARGS(1)
            success = fexpr_get_acb_raw(res, arg, prec);
            break;

        case FEXPR_Pow:
            ACB_FUNCTION_2(acb_pow)

        /* todo: acb_polygamma */
        case FEXPR_DigammaFunction:
            ACB_FUNCTION_1(acb_digamma)

        case FEXPR_Re:
            REQUIRE_NARGS(1)
            success = fexpr_get_acb_raw(res, arg, prec);
            arb_zero(acb_imagref(res));
            break;

        case FEXPR_RiemannZeta:
            ACB_FUNCTION_1(acb_zeta)

        case FEXPR_Root:
            ACB_FUNCTION_2(_acb_root)

        case FEXPR_RootOfUnity:
            if (nargs == 1 || nargs == 2)
            {
                fmpq_t q;
                fmpz_init(n);
                fmpz_init(m);
                fmpz_one(m);
                success = fexpr_get_fmpz(n, arg);
                success = success && (fmpz_sgn(n) >= 1);
                if (success && nargs == 2)
                {
                    fexpr_view_next(arg);
                    success = fexpr_get_fmpz(m, arg);
                }
                if (success)
                {
                    fmpz_mul_2exp(m, m, 1);
                    *fmpq_denref(q) = *n;
                    *fmpq_numref(q) = *m;
                    arb_sin_cos_pi_fmpq(acb_imagref(res), acb_realref(res), q, prec);
                }
            }
            break;

        case FEXPR_Sign:
            ACB_FUNCTION_1(acb_sgn)

        case FEXPR_Sin:
            ACB_FUNCTION_1(acb_sin)

        case FEXPR_Sinh:
            ACB_FUNCTION_1(acb_sinh)

        case FEXPR_Sqrt:
            ACB_FUNCTION_1(acb_sqrt)

        case FEXPR_Sub:
            if (nargs == 0)
            {
                acb_zero(res);
            }
            else if (nargs == 1)
            {
                success = fexpr_get_acb_raw(res, arg, prec);
            }
            else
            {
                acb_init(t);
                success = fexpr_get_acb_raw(res, arg, prec);
                for (i = 1; success && i < nargs; i++)
                {
                    fexpr_view_next(arg);
                    success = fexpr_get_acb_raw(t, arg, prec);
                    acb_sub(res, res, t, prec);
                }
                acb_clear(t);
            }
            break;

        case FEXPR_Tan:
            ACB_FUNCTION_1(acb_tan)

        case FEXPR_Tanh:
            ACB_FUNCTION_1(acb_tanh)

        default:
            success = 0;
        }

        if (!success || !acb_is_finite(res))
        {
            success = 0;
            acb_indeterminate(res);
        }

        return success;
    }
}

int
fexpr_get_acb_with_accuracy(acb_t res, const fexpr_t expr, slong prec, ulong flags)
{
    slong wp, initial, maxprec;
    int success = 0;

    initial = prec * 1.05 + 20;
    maxprec = FLINT_MAX(4096, 4 * initial);

    for (wp = initial; wp < maxprec; wp *= 2)
    {
        success = fexpr_get_acb_raw(res, expr, wp);

        if (acb_rel_accuracy_bits(res) >= prec)
            break;
    }

    return success;
}

char *
fexpr_get_decimal_str(const fexpr_t expr, slong digits, ulong flags)
{
    calcium_stream_t t;
    acb_t v;

    digits = FLINT_MAX(digits, 1);

    acb_init(v);
    calcium_stream_init_str(t);

    if (fexpr_get_acb_with_accuracy(v, expr, digits * 3.333 + 1, 0))
        calcium_write_acb(t, v, digits, ARB_STR_NO_RADIUS);
    else
        calcium_write(t, "?");

    acb_clear(v);
    return t->s;
}

/* todo: implement void fexpr_set_acb(fexpr_t res, const acb_t x); */
