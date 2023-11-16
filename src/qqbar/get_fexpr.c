/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "fmpq.h"
#include "fmpq_poly.h"
#include "arb_fmpz_poly.h"
#include "acb_poly.h"
#include "fexpr_builtin.h"
#include "ulong_extras.h"
#include "qqbar.h"

static ulong _deflation(const fmpz * poly, slong len)
{
    fmpz_poly_t t;
    t->alloc = t->length = len;
    t->coeffs = (fmpz *) poly;
    return arb_fmpz_poly_deflation(t);
}

void
fexpr_set_arb(fexpr_t res, const arb_t x)
{
    fexpr_t mid, rad, h;
    arf_t r;
    fexpr_init(mid);
    fexpr_init(rad);
    fexpr_init(h);

    fexpr_set_symbol_builtin(h, FEXPR_RealBall);
    fexpr_set_arf(mid, arb_midref(x));
    arf_init_set_mag_shallow(r, arb_radref(x));
    fexpr_set_arf(rad, r);

    fexpr_call2(res, h, mid, rad);

    fexpr_clear(mid);
    fexpr_clear(rad);
    fexpr_clear(h);
}

void
fexpr_set_acb(fexpr_t res, const acb_t x)
{
    if (arb_is_zero(acb_imagref(x)))
    {
        fexpr_set_arb(res, acb_realref(x));
    }
    else if (arb_is_zero(acb_realref(x)))
    {
        fexpr_t v, i;
        fexpr_init(v);
        fexpr_init(i);
        fexpr_set_arb(v, acb_imagref(x));
        fexpr_set_symbol_builtin(i, FEXPR_NumberI);
        fexpr_mul(res, v, i);
        fexpr_clear(v);
        fexpr_clear(i);
    }
    else
    {
        fexpr_t r, v, i;
        fexpr_init(r);
        fexpr_init(v);
        fexpr_init(i);
        fexpr_set_arb(v, acb_imagref(x));
        fexpr_set_symbol_builtin(i, FEXPR_NumberI);
        fexpr_mul(r, v, i);
        fexpr_set_arb(v, acb_realref(x));
        fexpr_add(res, v, r);
        fexpr_clear(r);
        fexpr_clear(v);
        fexpr_clear(i);
    }
}

/* todo: better code (no temporaries) */
void
fexpr_set_list_fmpz_poly(fexpr_t res, const fmpz_poly_t poly)
{
    fexpr_struct * coeffs;
    fexpr_t t;
    slong i, len;

    len = poly->length;

    coeffs = _fexpr_vec_init(len);
    fexpr_init(t);

    for (i = 0; i < len; i++)
        fexpr_set_fmpz(coeffs + i, poly->coeffs + i);

    fexpr_set_symbol_builtin(t, FEXPR_List);
    fexpr_call_vec(res, t, coeffs, len);

    _fexpr_vec_clear(coeffs, len);
    fexpr_clear(t);
}

void
qqbar_get_fexpr_repr(fexpr_t res, const qqbar_t x)
{
    slong j, d;
    /* todo: better fexpr construction code */
    fexpr_struct * coeffs;
    fexpr_t t, u, v, w;

    d = qqbar_degree(x);

    coeffs = _fexpr_vec_init(d + 1);
    fexpr_init(t);
    fexpr_init(u);
    fexpr_init(v);
    fexpr_init(w);

    for (j = 0; j <= d; j++)
        fexpr_set_fmpz(coeffs + j, QQBAR_COEFFS(x) + j);

    fexpr_set_symbol_builtin(u, FEXPR_List);
    fexpr_call_vec(t, u, coeffs, d + 1);
    fexpr_set_symbol_builtin(v, FEXPR_AlgebraicNumberSerialized);

    fexpr_set_acb(u, QQBAR_ENCLOSURE(x));

    fexpr_call2(res, v, t, u);

    _fexpr_vec_clear(coeffs, d + 1);
    fexpr_clear(t);
    fexpr_clear(u);
    fexpr_clear(v);
    fexpr_clear(w);
}


void
_qqbar_get_fexpr_root_nearest(fexpr_t res, const fmpz_poly_t poly, const char * re_s, const char * im_s)
{
    fexpr_t Decimal, a, b, I, s;

    fexpr_init(Decimal);
    fexpr_init(a);
    fexpr_init(b);
    fexpr_init(I);
    fexpr_init(s);

    fexpr_set_symbol_builtin(Decimal, FEXPR_Decimal);

    if (re_s == NULL && im_s == NULL)
    {
        fexpr_set_string(s, "0.0");
        fexpr_call1(a, Decimal, s);
    }
    else
    {
        if (re_s != NULL)
        {
            fexpr_set_string(s, re_s);
            fexpr_call1(a, Decimal, s);
        }

        if (im_s != NULL)
        {
            fexpr_set_string(s, im_s);
            fexpr_call1(b, Decimal, s);
            fexpr_set_symbol_builtin(I, FEXPR_NumberI);
            fexpr_mul(s, b, I);
            fexpr_swap(b, s);
        }
    }

    if (im_s == NULL)
        fexpr_swap(s, a);
    else if (re_s == NULL)
        fexpr_swap(s, b);
    else
        fexpr_add(s, a, b);

    fexpr_set_list_fmpz_poly(b, poly);
    fexpr_set_symbol_builtin(a, FEXPR_PolynomialRootNearest);

    fexpr_call2(res, a, b, s);

    fexpr_clear(Decimal);
    fexpr_clear(a);
    fexpr_clear(b);
    fexpr_clear(I);
    fexpr_clear(s);
}

void
qqbar_get_decimal_root_nearest(char ** re_s, char ** im_s, const qqbar_t x, slong default_digits)
{
    acb_t z, point, delta;
    acb_poly_t poly;
    arb_t lhs, rhs, R, Rpow, tmpr;
    slong k, d, digits, prec, wp;
    int success;
    int imag_zero, real_zero;

    d = qqbar_degree(x);

    /* For rational numbers, any numerical approximation will suffice */
    if (d == 1)
    {
        arb_t t;
        arb_init(t);
        digits = default_digits;
        qqbar_get_arb(t, x, 3.333 * digits + 10);
        re_s[0] = arb_get_str(t, digits, ARB_STR_NO_RADIUS);
        im_s[0] = NULL;
        arb_clear(t);
        return;
    }

    imag_zero = (qqbar_sgn_im(x) == 0);
    real_zero = (qqbar_sgn_re(x) == 0);

    /* For nonreal quadratics, getting the correct sign of the imaginary
       part (guaranteed by qqbar_get_acb) is sufficient. */
    if (d == 2 && !imag_zero)
    {
        acb_t t;
        acb_init(t);
        digits = default_digits;
        qqbar_get_acb(t, x, 3.333 * digits + 10);
        if (real_zero)
            re_s[0] = NULL;
        else
            re_s[0] = arb_get_str(acb_realref(t), digits, ARB_STR_NO_RADIUS);
        im_s[0] = arb_get_str(acb_imagref(t), digits, ARB_STR_NO_RADIUS);
        acb_clear(t);
        return;
    }

    acb_init(z);
    acb_init(point);
    acb_init(delta);
    acb_poly_init(poly);
    arb_init(lhs);
    arb_init(rhs);
    arb_init(R);
    arb_init(Rpow);
    arb_init(tmpr);

    acb_set(z, QQBAR_ENCLOSURE(x));
    if (imag_zero) arb_zero(acb_imagref(z));
    if (real_zero) arb_zero(acb_realref(z));

    re_s[0] = im_s[0] = NULL;

    success = 0;
    for (digits = default_digits; !success; digits *= 2)
    {
        prec = digits * 3.333 + 10;

        if (digits != 6 && 0)
            printf("digits " WORD_FMT "d\n", digits);

        if (acb_rel_accuracy_bits(z) < prec)
        {
            _qqbar_enclosure_raw(z, QQBAR_POLY(x), z, prec);
            if (imag_zero) arb_zero(acb_imagref(z));
            if (real_zero) arb_zero(acb_realref(z));
        }

        flint_free(re_s[0]);
        flint_free(im_s[0]);
        re_s[0] = arb_get_str(acb_realref(z), digits, ARB_STR_NO_RADIUS);
        im_s[0] = arb_get_str(acb_imagref(z), digits, ARB_STR_NO_RADIUS);

        /* Let x be the true root and let z = a + bi be the approximate point. */
        /* Verify that D(z, C*|z-x|) contains a unique root, for some C > 1. */

        for (wp = prec; ; wp *= 2)
        {
            /* printf("inner loop %ld / %ld,  acc %ld\n", wp, prec, acb_rel_accuracy_bits(z)); */

            if (acb_rel_accuracy_bits(z) < wp)
            {
                _qqbar_enclosure_raw(z, QQBAR_POLY(x), z, wp);
                if (imag_zero) arb_zero(acb_imagref(z));
                if (real_zero) arb_zero(acb_realref(z));
            }

            arb_set_str(acb_realref(point), re_s[0], wp);
            arb_set_str(acb_imagref(point), im_s[0], wp);

            acb_sub(delta, z, point, wp);
            acb_abs(R, delta, wp);

            /* need more precision */
            if (arb_contains_zero(R))
            {
                success = 0;
                continue;
            }

            /* C = 1.25 */
            arb_mul_ui(R, R, 5, wp);
            arb_mul_2exp_si(R, R, -2);

            acb_poly_set_fmpz_poly(poly, QQBAR_POLY(x), wp);
            acb_poly_taylor_shift(poly, poly, point, wp);

            /* |a_1| R > (pi/3) (|a_0| + |a_2| R^2 + ... + |a_d| R^d) */
            acb_abs(lhs, poly->coeffs + 1, wp);
            arb_mul(lhs, lhs, R, wp);

            acb_abs(rhs, poly->coeffs, wp);

            arb_set(Rpow, R);
            for (k = 2; k <= d; k++)
            {
                arb_mul(Rpow, Rpow, R, wp);
                acb_abs(tmpr, poly->coeffs + k, wp);
                arb_addmul(rhs, tmpr, Rpow, wp);
            }

            arb_const_pi(tmpr, wp);
            arb_mul(rhs, rhs, tmpr, wp);
            arb_div_ui(rhs, rhs, 3, wp);

            if (arb_overlaps(lhs, rhs))
                continue;

            success = arb_gt(lhs, rhs);
            break;
        }
    }

    if (real_zero)
    {
        flint_free(re_s[0]);
        re_s[0] = NULL;
    }

    if (imag_zero)
    {
        flint_free(im_s[0]);
        im_s[0] = NULL;
    }

    acb_clear(z);
    acb_clear(point);
    acb_clear(delta);
    acb_poly_clear(poly);
    arb_clear(lhs);
    arb_clear(rhs);
    arb_clear(R);
    arb_clear(Rpow);
    arb_clear(tmpr);
}

void
qqbar_get_fexpr_root_nearest(fexpr_t res, const qqbar_t x)
{
    char * re_s;
    char * im_s;

    qqbar_get_decimal_root_nearest(&re_s, &im_s, x, 6);
    _qqbar_get_fexpr_root_nearest(res, QQBAR_POLY(x), re_s, im_s);

    flint_free(re_s);
    flint_free(im_s);
}

void
qqbar_get_fexpr_root_indexed(fexpr_t res, const qqbar_t x)
{
    qqbar_ptr conjugates;
    slong d, i, j;

    d = qqbar_degree(x);

    conjugates = _qqbar_vec_init(d);
    qqbar_conjugates(conjugates, x);

    for (i = 0; i < d; i++)
    {
        if (qqbar_equal(conjugates + i, x))
        {
            /* todo: better fexpr construction code */
            fexpr_struct * coeffs;
            fexpr_t t, u, v;

            coeffs = _fexpr_vec_init(d + 1);
            fexpr_init(t);
            fexpr_init(u);
            fexpr_init(v);

            for (j = 0; j <= d; j++)
                fexpr_set_fmpz(coeffs + j, QQBAR_COEFFS(x) + j);

            fexpr_set_symbol_builtin(u, FEXPR_List);
            fexpr_call_vec(t, u, coeffs, d + 1);
            fexpr_set_si(u, i + 1);
            fexpr_set_symbol_builtin(v, FEXPR_PolynomialRootIndexed);
            fexpr_call2(res, v, t, u);

            _fexpr_vec_clear(coeffs, d + 1);
            fexpr_clear(t);
            fexpr_clear(u);
            fexpr_clear(v);

            break;
        }
    }

    _qqbar_vec_clear(conjugates, d);
}

void
fexpr_sqrt(fexpr_t res, const fexpr_t a)
{
    /* todo: handle aliasing in call1 */
    if (res == a)
    {
        fexpr_t tmp;
        fexpr_init(tmp);
        fexpr_set(tmp, a);
        fexpr_sqrt(res, tmp);
        fexpr_clear(tmp);
    }
    else
    {
        /* todo: avoid tmp alloc */
        fexpr_t tmp;
        fexpr_init(tmp);
        fexpr_set_symbol_builtin(tmp, FEXPR_Sqrt);
        fexpr_call1(res, tmp, a);
        fexpr_clear(tmp);
    }
}

void
fexpr_div_ui(fexpr_t res, const fexpr_t a, ulong c)
{
    fexpr_t t, u;
    fexpr_init(t);
    fexpr_init(u);
    fexpr_set_ui(u, c);
    fexpr_div(t, a, u);
    fexpr_swap(res, t);
    fexpr_clear(t);
    fexpr_clear(u);
}

/* cos(pi p/q) */
void
_fexpr_cos_pi_pq(fexpr_t res, slong p, ulong q)
{
    int sign = 1;
    int sine = 0;
    ulong g;
    fexpr_t t, u;

    if (p < 0)
    {
        _fexpr_cos_pi_pq(res, -p, q);
        return;
    }

    p = p % (2 * q);

    if (p > q)
    {
        p = 2 * q - p;
    }

    if (2 * p > q)
    {
        p = q - p;
        sign = -1;
    }

    if (p == 0)
    {
        fexpr_set_si(res, sign);
        return;
    }

    if (2 * p == q)
    {
        fexpr_set_ui(res, 0);
        return;
    }

    if (3 * p == q)
    {
        fexpr_set_si(res, sign);
        fexpr_div_ui(res, res, 2);
        return;
    }

    if (4 * p == q)
    {
        fexpr_set_ui(res, 2);
        fexpr_sqrt(res, res);
        fexpr_div_ui(res, res, 2);
        if (sign == -1)
            fexpr_neg(res, res);
        return;
    }

    if (6 * p == q)
    {
        fexpr_set_ui(res, 3);
        fexpr_sqrt(res, res);
        fexpr_div_ui(res, res, 2);
        if (sign == -1)
            fexpr_neg(res, res);
        return;
    }

    if (12 * p == q || 12 * p == 5 * q)
    {
        fexpr_init(t);
        fexpr_init(u);

        fexpr_set_ui(t, 3);
        fexpr_sqrt(t, t);
        fexpr_set_ui(u, 1);

        if (12 * p == q)
            fexpr_add(res, t, u);
        else
            fexpr_sub(res, t, u);

        fexpr_set_ui(t, 2);
        fexpr_sqrt(t, t);

        fexpr_mul(u, t, res);
        fexpr_div_ui(res, u, 4);

        if (sign == -1)
            fexpr_neg(res, res);

        fexpr_clear(t);
        fexpr_clear(u);
        return;
    }

    if (4 * p > q)
    {
        p = q - 2 * p;
        q = 2 * q;
        sine = 1;
    }

    g = n_gcd(p, q);
    if (g != 1)
    {
        p /= g;
        q /= g;
    }

    fexpr_init(t);
    fexpr_init(u);

    if (p == 1)
    {
        fexpr_set_symbol_builtin(res, FEXPR_Pi);
    }
    else
    {
        fexpr_set_ui(t, p);
        fexpr_set_symbol_builtin(u, FEXPR_Pi);
        fexpr_mul(res, t, u);
    }

    fexpr_div_ui(t, res, q);

    if (sine)
        fexpr_set_symbol_builtin(u, FEXPR_Sin);
    else
        fexpr_set_symbol_builtin(u, FEXPR_Cos);

    fexpr_call1(res, u, t);
    if (sign == -1)
        fexpr_neg(res, res);

    fexpr_clear(t);
    fexpr_clear(u);
}

ulong
qqbar_try_as_cyclotomic(qqbar_t zeta, fmpq_poly_t poly, const qqbar_t x)
{
    ulong * phi;
    ulong N1, N2, d2, order, d;
    slong p, q, i;
    double U;
    slong bits;
    int success;

    d = qqbar_degree(x);

    success = 0;
    order = 0;

    bits = 2 * qqbar_height_bits(x) /*+ 20 * d */ + 40;
    d2 = 4 * d;

    /* Compute inverse image of the totient function for d, 2*d and 4*d. */

    /* Determine lower and upper bounds [N1, N2) */
    U = d2;
    for (p = 2; p <= d2; p++)
        if (d2 % (p - 1) == 0 && n_is_prime(p))
            U = (U * p) / (p - 1);
    N1 = d + 1;
    N2 = U + 3;   /* +3 as safety for possible floating-point rounding */

/*        N2 = FLINT_MIN(N2, 1 + 30); */

    phi = flint_malloc(N2 * sizeof(ulong));

    for (i = 0; i < N2; i++)
        phi[i] = i;

    for (p = 2; p < N2; p++)
    {
        if (phi[p] == p)
        {
            phi[p] = p - 1;
            for (q = 2 * p; q < N2; q += p)
                phi[q] = (phi[q] / p) * (p - 1);
        }
    }

    for (i = N1; i < N2 && !success; i++)
    {
        if (phi[i] == d || phi[i] == 2 * d || phi[i] == 4 * d)
        {
            qqbar_root_of_unity(zeta, 1, i);

            if (qqbar_express_in_field(poly, zeta, x, bits, 0, bits))
            {
                order = i;
                success = 1;
                break;
            }
        }
    }

    flint_free(phi);

    return order;
}

/* poly(exp(2 pi i / n)) */
void
_qqbar_get_fexpr_cyclotomic(fexpr_t res, const fmpq_poly_t poly, slong n, int pure_real, int pure_imag)
{
    fexpr_vec_t terms;
    fexpr_t term, t, u, v, w;
    ulong p, q, g;
    slong i;

    fexpr_vec_init(terms, 0);
    fexpr_init(term);
    fexpr_init(t);
    fexpr_init(u);
    fexpr_init(v);
    fexpr_init(w);

    for (i = 0; i < poly->length; i++)
    {
        if (!fmpz_is_zero(poly->coeffs + i))
        {
            if (i == 0)
            {
                fexpr_set_fmpz(term, poly->coeffs + i);
            }
            else
            {
                p = 2 * i;
                q = n;
                g = n_gcd(p, q);
                p /= g;
                q /= g;

                if (pure_real)
                {
                    _fexpr_cos_pi_pq(v, p, q);
                }
                else
                {
                    fexpr_set_ui(t, p);
                    fexpr_set_symbol_builtin(u, FEXPR_Pi);
                    fexpr_set_symbol_builtin(v, FEXPR_NumberI);
                    fexpr_set_symbol_builtin(w, FEXPR_Mul);

                    if (p == 1)
                        fexpr_call2(term, w, u, v);
                    else
                        fexpr_call3(term, w, t, u, v);

                    fexpr_set_ui(t, q);
                    fexpr_div(u, term, t);

                    fexpr_set_symbol_builtin(w, FEXPR_Exp);
                    fexpr_call1(v, w, u);
                }

                if (fmpz_is_one(poly->coeffs + i))
                    fexpr_swap(term, v);
                else
                {
                    fexpr_set_fmpz(t, poly->coeffs + i);
                    fexpr_mul(term, t, v);
                }
            }

            fexpr_vec_append(terms, term);
        }
    }

    fexpr_set_symbol_builtin(t, FEXPR_Add);
    fexpr_call_vec(res, t, terms->entries, terms->length);

    if (!fmpz_is_one(poly->den))
    {
        fexpr_set_fmpz(t, poly->den);
        fexpr_div(u, res, t);
        fexpr_swap(res, u);
    }

    /* todo: also want this with expanded exponentials */
    if (pure_real)
    {
        fexpr_expanded_normal_form(res, res, 0);
    }

    fexpr_vec_clear(terms);
    fexpr_clear(term);
    fexpr_clear(t);
    fexpr_clear(u);
    fexpr_clear(v);
    fexpr_clear(w);
}

int
qqbar_get_fexpr_formula(fexpr_t res, const qqbar_t x, ulong flags)
{
    slong d;
    int success;

    d = qqbar_degree(x);

    if (d == 1)
    {
        fmpq_t r;
        fmpz_t t;
        fmpq_init(r);
        fmpz_init(t);
        qqbar_get_quadratic(fmpq_numref(r), t, t, fmpq_denref(r), x, 0);
        fexpr_set_fmpq(res, r);
        fmpq_clear(r);
        fmpz_clear(t);
        return 1;
    }

    if (d <= 2 && (flags & (QQBAR_FORMULA_QUADRATICS | QQBAR_FORMULA_GAUSSIANS)))
    {
        fmpz_t a, b, c, q;
        fexpr_t t, u, v;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        fmpz_init(q);
        fexpr_init(t);
        fexpr_init(u);
        fexpr_init(v);

        qqbar_get_quadratic(a, b, c, q, x, 2);

        if (fmpz_equal_si(c, -1) && (flags & QQBAR_FORMULA_GAUSSIANS))
        {
            fexpr_set_symbol_builtin(t, FEXPR_NumberI);
            success = 1;
        }
        else if (flags & QQBAR_FORMULA_QUADRATICS)
        {
            fexpr_set_fmpz(u, c);
            fexpr_set_symbol_builtin(v, FEXPR_Sqrt);
            fexpr_call1(t, v, u);
            success = 1;
        }
        else
        {
            success = 0;
        }

        if (success)
        {
            if (fmpz_is_zero(a))
            {
                if (fmpz_equal_si(b, -1))
                    fexpr_neg(u, t);
                else if (fmpz_equal_si(b, 1))
                    fexpr_swap(u, t);
                else
                {
                    fexpr_set_fmpz(v, b);
                    fexpr_mul(u, v, t);
                }
            }
            else
            {
                if (fmpz_equal_si(b, -1))
                {
                    fexpr_set_fmpz(v, a);
                    fexpr_sub(u, v, t);
                }
                else if (fmpz_equal_si(b, 1))
                {
                    fexpr_set_fmpz(v, a);
                    fexpr_add(u, v, t);
                }
                else
                {
                    fexpr_set_fmpz(u, b);
                    fexpr_mul(v, u, t);
                    fexpr_set_fmpz(t, a);
                    fexpr_add(u, t, v);
                }
            }

            if (fmpz_is_one(q))
            {
                fexpr_swap(res, u);
            }
            else
            {
                fexpr_set_fmpz(t, q);
                fexpr_div(res, u, t);
            }
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
        fmpz_clear(q);
        fexpr_clear(t);
        fexpr_clear(u);
        fexpr_clear(v);

        if (success)
            return 1;
    }

    /* Identify special constants */
    if (flags & QQBAR_FORMULA_CYCLOTOMICS)
    {
        slong p;
        ulong q;

        /* Quick detection of roots of unity. Todo: could also detect
           rational multiples of roots of unity. */
        if (qqbar_is_root_of_unity(&p, &q, x))
        {
            fexpr_t s, t, u, v, w;
            /* exp(2 p pi i / q) */

            if (q % 2 == 0)
                q /= 2;
            else
                p *= 2;

            fexpr_init(s);
            fexpr_init(t);
            fexpr_init(u);
            fexpr_init(v);
            fexpr_init(w);

            fexpr_set_symbol_builtin(s, FEXPR_Mul);
            fexpr_set_si(t, p);
            fexpr_set_symbol_builtin(u, FEXPR_Pi);
            fexpr_set_symbol_builtin(v, FEXPR_NumberI);
            if (p == 1)
                fexpr_call2(w, s, u, v);
            else
                fexpr_call3(w, s, t, u, v);

            fexpr_set_si(t, q);
            fexpr_div(u, w, t);
            fexpr_set_symbol_builtin(t, FEXPR_Exp);
            fexpr_call1(res, t, u);

            fexpr_clear(s);
            fexpr_clear(t);
            fexpr_clear(u);
            fexpr_clear(v);
            fexpr_clear(w);

            return 1;
        }
    }

    /* Check for elements of cyclotomic fields */
    if (flags & QQBAR_FORMULA_CYCLOTOMICS)
    {
        fmpq_poly_t poly;
        qqbar_t zeta;
        ulong q;

        fmpq_poly_init(poly);
        qqbar_init(zeta);

        q = qqbar_try_as_cyclotomic(zeta, poly, x);

        if (q != 0)
            _qqbar_get_fexpr_cyclotomic(res, poly, q, qqbar_sgn_im(x) == 0, qqbar_sgn_re(x) == 0);

        fmpq_poly_clear(poly);
        qqbar_clear(zeta);

        if (q != 0)
            return 1;
    }

    if (flags & QQBAR_FORMULA_DEFLATION)
    {
        slong deflation;
        ulong flags2;

        deflation = _deflation(QQBAR_COEFFS(x), d + 1);

        if (deflation > 1)
        {
            int neg;
            qqbar_t t, u, v;

            /* Avoid possible cases of infinite recursion. */
            /* Todo: in some cases, we might want to carefully allow deflation recursively. */
            flags2 = flags & ~QQBAR_FORMULA_DEFLATION;

            neg = qqbar_is_real(x) && qqbar_sgn_re(x) < 0;

            qqbar_init(t);
            qqbar_init(u);
            qqbar_init(v);

            if (neg)
                qqbar_neg(t, x);
            else
                qqbar_set(t, x);

            qqbar_pow_ui(u, t, deflation);
            success = qqbar_get_fexpr_formula(res, u, flags2);

            if (success)
            {
                fexpr_t a, b;

                fexpr_init(a);
                fexpr_init(b);

                /* todo: compute the root of unity more efficiently */
                qqbar_root_ui(v, u, deflation);
                qqbar_div(u, x, v);

                if (deflation == 2)
                {
                    fexpr_sqrt(b, res);
                }
                else
                {
                    fexpr_set_ui(b, 1);
                    fexpr_div_ui(a, b, deflation);
                    fexpr_pow(b, res, a);
                }

                if (qqbar_is_one(u))
                    fexpr_set(res, b);
                else if (qqbar_is_neg_one(u))
                    fexpr_neg(res, b);
                else
                {
                    /* Todo: | flags by GAUSSIANS, CYCLOTOMICS etc. to always express the root of unity? */
                    success = qqbar_get_fexpr_formula(a, u, flags2);

                    if (success)
                        fexpr_mul(res, b, a);
                }

                fexpr_clear(a);
                fexpr_clear(b);
            }

            qqbar_clear(t);
            qqbar_clear(u);
            qqbar_clear(v);

            if (success)
                return 1;
        }
    }

    /* Try a+bi or |x|*s separation. */
    if ((flags & QQBAR_FORMULA_SEPARATION) && !qqbar_is_real(x))
    {
        qqbar_t a, b;
        fexpr_t t, u, v;
        ulong flags2;

        qqbar_init(a);
        qqbar_init(b);
        fexpr_init(t);
        fexpr_init(u);
        fexpr_init(v);

        qqbar_re_im(a, b, x);

        /* Avoid possible cases of infinite recursion. */
        /* Todo: in some cases, we might want to carefully allow separation recursively,
           e.g. if the degree has been strictly decreased. */
        flags2 = flags & ~QQBAR_FORMULA_SEPARATION;

        success = (qqbar_degree(a) <= d) && (qqbar_degree(b) <= d) &&
                    qqbar_get_fexpr_formula(t, a, flags2) &&
                    qqbar_get_fexpr_formula(u, b, flags2);

        if (success)
        {
            fexpr_set_symbol_builtin(res, FEXPR_NumberI);
            fexpr_mul(v, u, res);

            if (qqbar_is_zero(a))
                fexpr_swap(res, v);
            else
                fexpr_add(res, t, v);
        }
        else
        {
            qqbar_abs(a, x);

            if (qqbar_degree(a) <= d && qqbar_get_fexpr_formula(t, a, flags2))
            {
                qqbar_div(b, x, a);

                if (qqbar_get_fexpr_formula(u, b, flags2))
                {
                    fexpr_mul(res, t, u);
                    success = 1;
                }
            }
        }

        qqbar_clear(a);
        qqbar_clear(b);
        fexpr_clear(t);
        fexpr_clear(u);
        fexpr_clear(v);

        if (success)
            return 1;
    }

    return 0;
}
