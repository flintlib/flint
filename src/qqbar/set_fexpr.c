/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "fmpz_poly.h"
#include "qqbar.h"
#include "fexpr.h"
#include "fexpr_builtin.h"

#ifdef __GNUC__
# define fabs __builtin_fabs
# define strchr __builtin_strchr
# define strlen __builtin_strlen
#else
# include <math.h>
# include <string.h>
#endif

int qqbar_set_fexpr(qqbar_t res, const fexpr_t expr);

int
_fexpr_parse_arf(arf_t res, const fexpr_t expr)
{
    if (fexpr_is_integer(expr))
    {
        fmpz_t m;
        fmpz_init(m);

        fexpr_get_fmpz(m, expr);
        arf_set_fmpz(res, m);

        fmpz_clear(m);
        return 1;
    }

    if (fexpr_is_builtin_call(expr, FEXPR_Neg) && fexpr_nargs(expr) == 1)
    {
        int success;
        fexpr_t t;
        fexpr_view_arg(t, expr, 0);
        success = _fexpr_parse_arf(res, t);
        arf_neg(res, res);
        return success;
    }

    if (fexpr_is_builtin_call(expr, FEXPR_Div) && fexpr_nargs(expr) == 2)
    {
        int success;
        fexpr_t num, den;
        fmpz_t p, q;

        fexpr_view_arg(num, expr, 0);
        fexpr_view_arg(den, expr, 1);

        fmpz_init(p);
        fmpz_init(q);

        success = fexpr_get_fmpz(p, num) && fexpr_get_fmpz(q, den);

        success = success && (fmpz_sgn(q) > 0) && (fmpz_val2(q) == fmpz_bits(q) - 1);

        if (success)
        {
            arf_set_fmpz(res, p);
            arf_mul_2exp_si(res, res, -(fmpz_bits(q) - 1));
        }

        fmpz_clear(p);
        fmpz_clear(q);

        return success;
    }

    if (fexpr_is_builtin_call(expr, FEXPR_Pow) && fexpr_nargs(expr) == 2)
    {
        fexpr_t base, exp;

        fexpr_view_arg(base, expr, 0);
        fexpr_view_arg(exp, expr, 1);

        if (fexpr_equal_ui(base, 2))
        {
            fmpz_t m, e;
            int success;

            fmpz_init(m);
            fmpz_init(e);

            success = fexpr_get_fmpz(e, exp);

            if (success)
            {
                arf_one(res);
                arf_mul_2exp_fmpz(res, res, e);
            }

            fmpz_clear(m);
            fmpz_clear(e);

            return success;
        }
    }

    if (fexpr_is_builtin_call(expr, FEXPR_Mul) && fexpr_nargs(expr) == 2)
    {
        fexpr_t man, pow, base, exp;

        fexpr_view_arg(man, expr, 0);
        fexpr_view_arg(pow, expr, 1);

        if (fexpr_is_builtin_call(pow, FEXPR_Pow) && fexpr_nargs(expr) == 2)
        {
            fexpr_view_arg(base, pow, 0);
            fexpr_view_arg(exp, pow, 1);

            if (fexpr_equal_ui(base, 2))
            {
                fmpz_t m, e;
                int success;

                fmpz_init(m);
                fmpz_init(e);

                success = fexpr_get_fmpz(m, man) && fexpr_get_fmpz(e, exp);

                if (success)
                {
                    arf_set_fmpz(res, m);
                    arf_mul_2exp_fmpz(res, res, e);
                }

                fmpz_clear(m);
                fmpz_clear(e);

                return success;
            }
        }
    }

    return 0;
}

int
_fexpr_parse_mag(mag_t res, const fexpr_t expr)
{
    int success;
    arf_t t;
    arf_init(t);
    success = _fexpr_parse_arf(t, expr);

    success = success && (arf_sgn(t) >= 0) && arf_is_finite(t) && (arf_bits(t) <= MAG_BITS);

    if (success)
    {
        fmpz_t m, e;

        fmpz_init(m);
        fmpz_init(e);

        arf_get_fmpz_2exp(m, e, t);

        mag_set_ui(res, fmpz_get_ui(m));
        mag_mul_2exp_fmpz(res, res, e);

        fmpz_clear(m);
        fmpz_clear(e);
    }

    arf_clear(t);
    return success;
}


int
_fexpr_parse_arb(arb_t res, const fexpr_t expr)
{
    if (fexpr_is_builtin_call(expr, FEXPR_RealBall) && fexpr_nargs(expr) == 2)
    {
        fexpr_t t, u;

        fexpr_view_arg(t, expr, 0);
        fexpr_view_arg(u, expr, 1);

        return _fexpr_parse_arf(arb_midref(res), t) && _fexpr_parse_mag(arb_radref(res), u);
    }

    return 0;
}

int
_fexpr_parse_acb(acb_t res, const fexpr_t expr)
{
    fexpr_t t, u;

    if (fexpr_is_builtin_call(expr, FEXPR_RealBall) && fexpr_nargs(expr) == 2)
    {
        arb_zero(acb_imagref(res));
        return _fexpr_parse_arb(acb_realref(res), expr);
    }

    if (fexpr_is_builtin_call(expr, FEXPR_Mul) && fexpr_nargs(expr) == 2)
    {
        fexpr_view_arg(t, expr, 1);
        if (!fexpr_is_builtin_symbol(t, FEXPR_NumberI))
            return 0;

        fexpr_view_arg(u, expr, 0);
        arb_zero(acb_realref(res));
        return _fexpr_parse_arb(acb_imagref(res), u);
    }

    if (fexpr_is_builtin_call(expr, FEXPR_Add) && fexpr_nargs(expr) == 2)
    {
        fexpr_view_arg(t, expr, 0);
        fexpr_view_arg(u, expr, 1);

        if (_fexpr_parse_acb(res, u) && arb_is_zero(acb_realref(res)))
            return _fexpr_parse_arb(acb_realref(res), t);
    }

    return 0;
}

int
fmpq_set_decimal(fmpq_t res, const char * inp, slong max_bits)
{
    char * emarker;
    char * buf;
    int success;
    slong i;
    fmpz_t exp;
    fmpz_t man;
    slong num_int, num_frac;
    int after_radix;

    if (inp[0] == '+')
    {
        return fmpq_set_decimal(res, inp + 1, max_bits);
    }

    if (inp[0] == '-')
    {
        success = fmpq_set_decimal(res, inp + 1, max_bits);
        fmpq_neg(res, res);
        return success;
    }

    success = 1;
    fmpz_init(exp);
    fmpz_init(man);
    buf = flint_malloc(strlen(inp) + 1);

    emarker = strchr(inp, 'e');
    if (emarker == NULL)
        emarker = strchr(inp, 'E');

    /* parse exponent (0 by default) */
    if (emarker != NULL)
    {
        /* allow e+42 as well as e42 */
        if (emarker[1] == '+')
        {
            if (!(emarker[2] >= '0' && emarker[2] <= '9'))
                success = 0;
            else
                success = !fmpz_set_str(exp, emarker + 2, 10);
        }
        else
            success = !fmpz_set_str(exp, emarker + 1, 10);

        if (!success)
            goto cleanup;
    }

    /* parse floating-point part */
    {
        num_int = 0;
        num_frac = 0;
        after_radix = 0;

        for (i = 0; inp + i != emarker && inp[i] != '\0'; i++)
        {
            if (inp[i] == '.' && !after_radix)
            {
                after_radix = 1;
            }
            else if (inp[i] >= '0' && inp[i] <= '9')
            {
                buf[num_int + num_frac] = inp[i];

                num_frac += after_radix;
                num_int += !after_radix;
            }
            else
            {
                success = 0;
                goto cleanup;
            }
        }

        buf[num_int + num_frac] = '\0';

        /* put trailing zeros into the exponent */
        while (num_int + num_frac > 1 && buf[num_int + num_frac - 1] == '0')
        {
            buf[num_int + num_frac - 1] = '\0';
            num_frac--;
        }

        fmpz_sub_si(exp, exp, num_frac);

        success = !fmpz_set_str(man, buf, 10);
        if (!success)
            goto cleanup;
    }

    if (fmpz_is_zero(man))
    {
        fmpq_zero(res);
    }
    else if (COEFF_IS_MPZ(*exp))
    {
        success = 0;
    }
    else
    {
        slong e = *exp;
        double size;

        size = fmpz_bits(man) + e * 3.321928094887;
        size = fabs(size);

        if (size > max_bits)
        {
            success = 0;
        }
        else
        {
            if (e >= 0)
            {
                fmpz_set_ui(exp, 10);
                fmpz_pow_ui(exp, exp, e);
                fmpz_mul(fmpq_numref(res), man, exp);
                fmpz_one(fmpq_denref(res));
            }
            else
            {
                fmpz_set_ui(exp, 10);
                fmpz_pow_ui(exp, exp, -e);
                fmpz_set(fmpq_numref(res), man);
                fmpz_set(fmpq_denref(res), exp);
                fmpq_canonicalise(res);
            }
        }
    }

    if (!success)
    {
        fmpq_zero(res);
    }

cleanup:
    fmpz_clear(exp);
    fmpz_clear(man);
    flint_free(buf);

    return success;
}

static int
_fexpr_check_pi_in_product(const fexpr_t expr)
{
    fexpr_t func, arg, arg2;
    slong i, nargs;
    int status, arg_status;

    if (fexpr_is_atom(expr))
    {
        if (fexpr_is_builtin_symbol(expr, FEXPR_Pi))
            return 1;

        return 0;
    }

    nargs = fexpr_nargs(expr);
    fexpr_view_func(func, expr);

    if (nargs == 1 && (fexpr_is_builtin_symbol(func, FEXPR_Neg) || fexpr_is_builtin_symbol(func, FEXPR_Pos)))
    {
        fexpr_view_arg(arg, expr, 0);
        return _fexpr_check_pi_in_product(arg);
    }

    if (nargs == 2 && (fexpr_is_builtin_symbol(func, FEXPR_Div)))
    {
        fexpr_view_arg(arg, expr, 0);
        fexpr_view_arg(arg2, expr, 1);
        if (_fexpr_check_pi_in_product(arg2) != 0)
            return -1;
        return _fexpr_check_pi_in_product(arg);
    }

    if (nargs >= 1 && (fexpr_is_builtin_symbol(func, FEXPR_Mul)))
    {
        status = 0;
        fexpr_view_arg(arg, expr, 0);
        for (i = 0; i < nargs; i++)
        {
            arg_status = _fexpr_check_pi_in_product(arg);
            if (arg_status == -1)
                return -1;
            if (arg_status == 1 && status == 1)
                return -1;
            if (arg_status == 1)
                status = 1;
            fexpr_view_next(arg);
        }

        return status;
    }

    return -1;
}

static int
_fexpr_get_rational_arg_pi(fmpq_t res, const fexpr_t expr, int times_i)
{
    int status, success;

    status = _fexpr_check_pi_in_product(expr);

    if (status == 0 || status == 1)
    {
        fexpr_t tmp, pi, one;
        qqbar_t v, i;
        fexpr_init(tmp);
        fexpr_init(pi);
        fexpr_init(one);
        qqbar_init(v);

        fexpr_set_symbol_builtin(pi, FEXPR_Pi);
        fexpr_set_si(one, 1);
        fexpr_replace(tmp, expr, pi, one);
        success = qqbar_set_fexpr(v, tmp);
        if (success)
        {
            if (times_i)
            {
                qqbar_init(i);
                qqbar_i(i);
                qqbar_div(v, v, i);
                qqbar_clear(i);
            }

            success = qqbar_is_rational(v);

            if (success)
            {
                fmpz_neg(fmpq_numref(res), QQBAR_COEFFS(v));
                fmpz_set(fmpq_denref(res), QQBAR_COEFFS(v) + 1);
            }
        }

        fexpr_clear(tmp);
        fexpr_clear(pi);
        fexpr_clear(one);
        qqbar_clear(v);

        return success;
    }

    return 0;
}

void
qqbar_set_fmpz_poly_root_indexed(qqbar_t res, const fmpz_poly_t poly, slong root_index)
{
    qqbar_ptr roots;
    slong d;
    d = fmpz_poly_degree(poly);
    roots = _qqbar_vec_init(d);
    qqbar_roots_fmpz_poly(roots, poly, 0);
    qqbar_set(res, roots + root_index - 1);
    _qqbar_vec_clear(roots, d);
}

void
qqbar_set_fmpz_poly_root_nearest(qqbar_t res, const fmpz_poly_t poly, const qqbar_t point)
{
    qqbar_ptr roots;
    slong i, best, d;
    acb_t t;
    arb_t distance, best_distance;
    int overlapping;

    d = fmpz_poly_degree(poly);

    roots = _qqbar_vec_init(d);
    acb_init(t);
    arb_init(distance);
    arb_init(best_distance);

    qqbar_roots_fmpz_poly(roots, poly, 0);

    acb_sub(t, QQBAR_ENCLOSURE(point), QQBAR_ENCLOSURE(roots), QQBAR_DEFAULT_PREC);
    acb_abs(best_distance, t, QQBAR_DEFAULT_PREC);
    best = 0;
    overlapping = 0;

    for (i = 1; i < d; i++)
    {
        acb_sub(t, QQBAR_ENCLOSURE(point), QQBAR_ENCLOSURE(roots + i), QQBAR_DEFAULT_PREC);
        acb_abs(distance, t, QQBAR_DEFAULT_PREC);

        if (arb_lt(distance, best_distance))
        {
            arb_swap(best_distance, distance);
            best = i;
            overlapping = 0;
        }
        else if (arb_overlaps(distance, best_distance))
        {
            overlapping = 1;
        }
    }

    if (overlapping)
    {
        qqbar_t exact_distance, best_exact_distance;

        qqbar_init(exact_distance);
        qqbar_init(best_exact_distance);

        qqbar_sub(best_exact_distance, point, roots + best);
        qqbar_abs2(best_exact_distance, best_exact_distance);

        for (i = 0; i < d; i++)
        {
            if (i != best)
            {
                acb_sub(t, QQBAR_ENCLOSURE(point), QQBAR_ENCLOSURE(roots + i), QQBAR_DEFAULT_PREC);
                acb_abs(distance, t, QQBAR_DEFAULT_PREC);

                if (arb_gt(distance, best_distance))
                    continue;

                qqbar_sub(exact_distance, point, roots + i);
                qqbar_abs2(exact_distance, exact_distance);

                if (qqbar_cmp_re(exact_distance, best_exact_distance) < 0)
                {
                    qqbar_swap(best_exact_distance, exact_distance);
                    best = i;
                }
            }
        }

        qqbar_clear(exact_distance);
        qqbar_clear(best_exact_distance);
    }

    qqbar_swap(res, roots + best);

    acb_clear(t);
    arb_clear(distance);
    arb_clear(best_distance);
    _qqbar_vec_clear(roots, d);
}


int
qqbar_set_fexpr(qqbar_t res, const fexpr_t expr)
{
    fexpr_t func, arg;
    slong id, i, nargs;
    qqbar_t t;
    fmpq_t q;
    int success;

    if (fexpr_is_integer(expr))
    {
        fmpz_t t;
        fmpz_init(t);
        fexpr_get_fmpz(t, expr);
        qqbar_set_fmpz(res, t);
        fmpz_clear(t);
        return 1;
    }

    if (fexpr_is_atom(expr))
    {
        if (fexpr_is_builtin_symbol(expr, FEXPR_NumberI))
        {
            qqbar_i(res);
            return 1;
        }

        if (fexpr_is_builtin_symbol(expr, FEXPR_GoldenRatio))
        {
            qqbar_phi(res);
            return 1;
        }

        return 0;
    }

    nargs = fexpr_nargs(expr);
    fexpr_view_func(func, expr);

    id = FEXPR_BUILTIN_ID(func->data[0]);

    switch (id)
    {
        case FEXPR_AlgebraicNumberSerialized:
        case FEXPR_PolynomialRootIndexed:
        case FEXPR_PolynomialRootNearest:
            if (nargs == 2)
            {
                slong root_index = 0;

                fexpr_view_arg(arg, expr, 1);

                if (id == FEXPR_PolynomialRootIndexed)
                {
                    fmpz_t tmp;
                    fmpz_init(tmp);
                    success = fexpr_get_fmpz(tmp, arg);
                    root_index = *tmp;
                    success = success && (root_index >= 1 && root_index <= COEFF_MAX);
                    fmpz_clear(tmp);
                }
                else if (id == FEXPR_PolynomialRootNearest)
                {
                    root_index = 1;
                    success = qqbar_set_fexpr(res, arg);
                }
                else
                {
                    success = _fexpr_parse_acb(QQBAR_ENCLOSURE(res), arg);
                }

                fexpr_view_arg(arg, expr, 0);
                success = success && fexpr_is_builtin_call(arg, FEXPR_List);

                if (success)
                {
                    slong deg, len;
                    len = fexpr_nargs(arg);
                    deg = len - 1;
                    success = (deg >= 1) && (deg >= root_index);

                    if (success)
                    {
                        fmpz_poly_t poly;
                        fexpr_t c;

                        fmpz_poly_init2(poly, len);
                        _fmpz_poly_set_length(poly, len);
                        fexpr_view_arg(c, arg, 0);
                        for (i = 0; i < len && success; i++)
                        {
                            success = fexpr_get_fmpz(poly->coeffs + i, c);
                            fexpr_view_next(c);
                        }

                        if (success)
                        {
                            _fmpz_poly_normalise(poly);
                            deg = fmpz_poly_degree(poly);

                            success = (deg >= root_index);

                            if (success && id == FEXPR_PolynomialRootIndexed)
                                qqbar_set_fmpz_poly_root_indexed(res, poly, root_index);
                            else if (success && id == FEXPR_PolynomialRootNearest)
                                qqbar_set_fmpz_poly_root_nearest(res, poly, res);
                            else if (success)
                                fmpz_poly_swap(QQBAR_POLY(res), poly);
                        }

                        fmpz_poly_clear(poly);

                        return success;
                    }
                }
            }
            break;

        case FEXPR_Decimal:
            if (nargs == 1)
            {
                fexpr_view_arg(arg, expr, 0);
                success = fexpr_is_string(arg);

                if (success)
                {
                    char * s = fexpr_get_string(arg);
                    fmpq_init(q);
                    success = fmpq_set_decimal(q, s, COEFF_MAX);
                    if (success)
                    {
                        qqbar_set_fmpq(res, q);
                    }
                    flint_free(s);
                    fmpq_clear(q);
                }

                return success;
            }
            break;

        case FEXPR_RootOfUnity:
            if (nargs == 1 || nargs == 2)
            {
                fmpz_t n, k;
                fmpz_init(n);
                fmpz_init(k);
                qqbar_init(t);
                fmpz_one(k);

                fexpr_view_arg(arg, expr, 0);
                success = qqbar_set_fexpr(t, arg) && qqbar_is_integer(t) && (qqbar_sgn_re(t) == 1);

                if (success)
                    fmpz_neg(n, QQBAR_COEFFS(t));

                success = success && !COEFF_IS_MPZ(*n);

                if (success && nargs == 2)
                {
                    fexpr_view_arg(arg, expr, 1);
                    success = qqbar_set_fexpr(res, arg) && qqbar_is_integer(res);
                    if (success)
                    {
                        fmpz_neg(k, QQBAR_COEFFS(res));
                        fmpz_fdiv_r(k, k, n);
                    }
                }

                if (success)
                    qqbar_root_of_unity(res, *k, *n);

                fmpz_clear(n);
                fmpz_clear(k);
                qqbar_clear(t);
                return success;
            }
            break;

        case FEXPR_Exp:
            if (nargs == 1)
            {
                fexpr_view_arg(arg, expr, 0);
                fmpq_init(q);
                success = _fexpr_get_rational_arg_pi(q, arg, 1);
                if (success && !COEFF_IS_MPZ(*fmpq_numref(q)) && !COEFF_IS_MPZ(*fmpq_denref(q)))
                    qqbar_exp_pi_i(res, *fmpq_numref(q), *fmpq_denref(q));
                fmpq_clear(q);
                return success;
            }
            break;

        case FEXPR_Sin:
            if (nargs == 1)
            {
                fexpr_view_arg(arg, expr, 0);
                fmpq_init(q);
                success = _fexpr_get_rational_arg_pi(q, arg, 0);
                if (success && !COEFF_IS_MPZ(*fmpq_numref(q)) && !COEFF_IS_MPZ(*fmpq_denref(q)))
                    qqbar_sin_pi(res, *fmpq_numref(q), *fmpq_denref(q));
                fmpq_clear(q);
                return success;
            }
            break;

        case FEXPR_Cos:
            if (nargs == 1)
            {
                fexpr_view_arg(arg, expr, 0);
                fmpq_init(q);
                success = _fexpr_get_rational_arg_pi(q, arg, 0);
                if (success && !COEFF_IS_MPZ(*fmpq_numref(q)) && !COEFF_IS_MPZ(*fmpq_denref(q)))
                    qqbar_cos_pi(res, *fmpq_numref(q), *fmpq_denref(q));
                fmpq_clear(q);
                return success;
            }
            break;

        case FEXPR_Tan:
            if (nargs == 1)
            {
                fexpr_view_arg(arg, expr, 0);
                fmpq_init(q);
                success = _fexpr_get_rational_arg_pi(q, arg, 0);
                if (success && !COEFF_IS_MPZ(*fmpq_numref(q)) && !COEFF_IS_MPZ(*fmpq_denref(q)))
                    success = qqbar_tan_pi(res, *fmpq_numref(q), *fmpq_denref(q));
                fmpq_clear(q);
                return success;
            }
            break;

        case FEXPR_Cot:
            if (nargs == 1)
            {
                fexpr_view_arg(arg, expr, 0);
                fmpq_init(q);
                success = _fexpr_get_rational_arg_pi(q, arg, 0);
                if (success && !COEFF_IS_MPZ(*fmpq_numref(q)) && !COEFF_IS_MPZ(*fmpq_denref(q)))
                    success = qqbar_cot_pi(res, *fmpq_numref(q), *fmpq_denref(q));
                fmpq_clear(q);
                return success;
            }
            break;

        case FEXPR_Sec:
            if (nargs == 1)
            {
                fexpr_view_arg(arg, expr, 0);
                fmpq_init(q);
                success = _fexpr_get_rational_arg_pi(q, arg, 0);
                if (success && !COEFF_IS_MPZ(*fmpq_numref(q)) && !COEFF_IS_MPZ(*fmpq_denref(q)))
                    success = qqbar_sec_pi(res, *fmpq_numref(q), *fmpq_denref(q));
                fmpq_clear(q);
                return success;
            }
            break;

        case FEXPR_Csc:
            if (nargs == 1)
            {
                fexpr_view_arg(arg, expr, 0);
                fmpq_init(q);
                success = _fexpr_get_rational_arg_pi(q, arg, 0);
                if (success && !COEFF_IS_MPZ(*fmpq_numref(q)) && !COEFF_IS_MPZ(*fmpq_denref(q)))
                    success = qqbar_csc_pi(res, *fmpq_numref(q), *fmpq_denref(q));
                fmpq_clear(q);
                return success;
            }
            break;

        case FEXPR_Pos:
            if (nargs == 1)
            {
                fexpr_view_arg(arg, expr, 0);
                success = qqbar_set_fexpr(res, arg);
                return success;
            }
            break;

        case FEXPR_Neg:
            if (nargs == 1)
            {
                fexpr_view_arg(arg, expr, 0);
                success = qqbar_set_fexpr(res, arg);
                if (success)
                    qqbar_neg(res, res);
                return success;
            }
            break;

        case FEXPR_Sqrt:
            if (nargs == 1)
            {
                fexpr_view_arg(arg, expr, 0);
                success = qqbar_set_fexpr(res, arg);
                if (success)
                    qqbar_sqrt(res, res);
                return success;
            }
            break;

        case FEXPR_Conjugate:
            if (nargs == 1)
            {
                fexpr_view_arg(arg, expr, 0);
                success = qqbar_set_fexpr(res, arg);
                if (success)
                    qqbar_conj(res, res);
                return success;
            }
            break;

        case FEXPR_Re:
            if (nargs == 1)
            {
                fexpr_view_arg(arg, expr, 0);
                success = qqbar_set_fexpr(res, arg);
                if (success)
                    qqbar_re(res, res);
                return success;
            }
            break;

        case FEXPR_Im:
            if (nargs == 1)
            {
                fexpr_view_arg(arg, expr, 0);
                success = qqbar_set_fexpr(res, arg);
                if (success)
                    qqbar_im(res, res);
                return success;
            }
            break;

        case FEXPR_Floor:
            if (nargs == 1)
            {
                fexpr_view_arg(arg, expr, 0);
                success = qqbar_set_fexpr(res, arg);
                if (success)
                {
                    fmpz_t n;
                    fmpz_init(n);
                    qqbar_floor(n, res);
                    qqbar_set_fmpz(res, n);
                    fmpz_clear(n);
                }
                return success;
            }
            break;

        case FEXPR_Ceil:
            if (nargs == 1)
            {
                fexpr_view_arg(arg, expr, 0);
                success = qqbar_set_fexpr(res, arg);
                if (success)
                {
                    fmpz_t n;
                    fmpz_init(n);
                    qqbar_ceil(n, res);
                    qqbar_set_fmpz(res, n);
                    fmpz_clear(n);
                }
                return success;
            }
            break;

        case FEXPR_Abs:
            if (nargs == 1)
            {
                fexpr_view_arg(arg, expr, 0);
                success = qqbar_set_fexpr(res, arg);
                if (success)
                    qqbar_abs(res, res);
                return success;
            }
            break;

        case FEXPR_Sign:
            if (nargs == 1)
            {
                fexpr_view_arg(arg, expr, 0);
                success = qqbar_set_fexpr(res, arg);
                if (success)
                    qqbar_sgn(res, res);
                return success;
            }
            break;

        case FEXPR_Csgn:
            if (nargs == 1)
            {
                fexpr_view_arg(arg, expr, 0);
                success = qqbar_set_fexpr(res, arg);
                if (success)
                    qqbar_set_si(res, qqbar_csgn(res));
                return success;
            }
            break;

        case FEXPR_Sub:
            if (nargs == 2)
            {
                qqbar_init(t);
                fexpr_view_arg(arg, expr, 0);
                success = qqbar_set_fexpr(res, arg);
                if (success)
                {
                    fexpr_view_next(arg);
                    success = qqbar_set_fexpr(t, arg);
                    if (success)
                        qqbar_sub(res, res, t);
                }
                qqbar_clear(t);
                return success;
            }
            break;

        case FEXPR_Div:
            if (nargs == 2)
            {
                qqbar_init(t);
                fexpr_view_arg(arg, expr, 0);
                success = qqbar_set_fexpr(res, arg);
                if (success)
                {
                    fexpr_view_next(arg);
                    success = qqbar_set_fexpr(t, arg);
                    success = success && !qqbar_is_zero(t);
                    if (success)
                        qqbar_div(res, res, t);
                }
                qqbar_clear(t);
                return success;
            }
            break;

        case FEXPR_Pow:
            if (nargs == 2)
            {
                qqbar_init(t);
                fexpr_view_arg(arg, expr, 0);
                success = qqbar_set_fexpr(res, arg);
                if (success)
                {
                    fexpr_view_next(arg);
                    success = qqbar_set_fexpr(t, arg);
                    success = success && qqbar_is_rational(t) && (!qqbar_is_zero(res) || qqbar_sgn_re(t) >= 0);
                    if (success)
                    {
                        fmpz_t p, q;
                        fmpz_init(p);
                        fmpz_init(q);
                        fmpz_neg(p, QQBAR_COEFFS(t));
                        fmpz_set(q, QQBAR_COEFFS(t) + 1);

                        success = (fmpz_bits(q) <= 20 && fmpz_bits(p) <= FLINT_BITS - 4);

                        if (success)
                        {
                            qqbar_root_ui(res, res, *q);
                            if (*p >= 0)
                                qqbar_pow_ui(res, res, *p);
                            else
                            {
                                qqbar_pow_ui(res, res, -(*p));
                                qqbar_inv(res, res);
                            }
                        }
                    }
                }
                qqbar_clear(t);
                return success;
            }
            break;

        case FEXPR_Add:
            if (nargs == 0)
            {
                qqbar_zero(res);
                return 1;
            }

            fexpr_view_arg(arg, expr, 0);
            success = qqbar_set_fexpr(res, arg);

            if (success && nargs > 1)
            {
                qqbar_init(t);
                for (i = 1; i < nargs && success; i++)
                {
                    fexpr_view_next(arg);
                    success = qqbar_set_fexpr(t, arg);
                    if (success)
                        qqbar_add(res, res, t);
                }
                qqbar_clear(t);
            }

            return success;

        case FEXPR_Mul:
            if (nargs == 0)
            {
                qqbar_one(res);
                return 1;
            }

            fexpr_view_arg(arg, expr, 0);
            success = qqbar_set_fexpr(res, arg);

            if (success && nargs > 1)
            {
                qqbar_init(t);
                for (i = 1; i < nargs && success; i++)
                {
                    fexpr_view_next(arg);
                    success = qqbar_set_fexpr(t, arg);
                    if (success)
                        qqbar_mul(res, res, t);
                }
                qqbar_clear(t);
            }

            return success;

        case FEXPR_Min:
        case FEXPR_Max:
            if (nargs >= 1)
            {
                fexpr_view_arg(arg, expr, 0);
                success = qqbar_set_fexpr(res, arg) && qqbar_is_real(res);

                if (success && nargs >= 2)
                {
                    qqbar_init(t);
                    for (i = 1; i < nargs && success; i++)
                    {
                        fexpr_view_next(arg);
                        success = qqbar_set_fexpr(t, arg) && qqbar_is_real(t);
                        if (success)
                        {
                            if ((qqbar_cmp_re(res, t) < 0) == (id == FEXPR_Max))
                                qqbar_swap(res, t);
                        }
                    }
                    qqbar_clear(t);
                }

                return success;
            }

        default:
            break;
    }

    return 0;
}
