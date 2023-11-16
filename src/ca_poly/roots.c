/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_poly.h"

void
_ca_poly_roots_quadratic(ca_t r1, ca_t r2, const ca_t a, const ca_t b, const ca_t c, ca_ctx_t ctx)
{
    ca_t d, t;

    ca_init(d, ctx);
    ca_init(t, ctx);

    ca_mul(t, a, c, ctx);
    ca_mul_ui(t, t, 4, ctx);
    ca_sqr(d, b, ctx);
    ca_sub(d, d, t, ctx);
    ca_sqrt(d, d, ctx);

    ca_inv(t, a, ctx);
    ca_div_ui(t, t, 2, ctx);

    ca_sub(r1, d, b, ctx);
    ca_add(r2, b, d, ctx);
    ca_neg(r2, r2, ctx);
    ca_mul(r1, r1, t, ctx);
    ca_mul(r2, r2, t, ctx);

    ca_clear(d, ctx);
    ca_clear(t, ctx);
}

/* exp(2 pi i * (1 / 3)) */
/* todo: make this faster to construct */
void
ca_omega(ca_t res, ca_ctx_t ctx)
{
    ca_pi_i(res, ctx);
    ca_mul_ui(res, res, 2, ctx);
    ca_div_ui(res, res, 3, ctx);
    ca_exp(res, res, ctx);
}

/* Solves a cubic equation using the cubic formula.
   Assumes (does not check) that a is invertible. */
int
_ca_poly_roots_cubic(ca_t r1, ca_t r2, ca_t r3, const ca_t a, const ca_t b, const ca_t c, const ca_t d, ca_ctx_t ctx)
{
    ca_t D0, D1, C, w1, w2, t;
    int success;

    ca_init(D0, ctx);
    ca_init(D1, ctx);
    ca_init(C, ctx);
    ca_init(w1, ctx);
    ca_init(w2, ctx);
    ca_init(t, ctx);

    /* D0 = b^2 - 3*a*c */
    ca_sqr(D0, b, ctx);
    ca_mul(t, a, c, ctx);
    ca_mul_ui(t, t, 3, ctx);
    ca_sub(D0, D0, t, ctx);

    /* D1 = b*(2*b^2 - 9*a*c) + 27*a^2*d */
    ca_sqr(D1, b, ctx);
    ca_mul_ui(D1, D1, 2, ctx);
    ca_mul(t, a, c, ctx);
    ca_mul_ui(t, t, 9, ctx);
    ca_sub(D1, D1, t, ctx);
    ca_mul(D1, D1, b, ctx);
    ca_sqr(t, a, ctx);
    ca_mul(t, t, d, ctx);
    ca_mul_ui(t, t, 27, ctx);
    ca_add(D1, D1, t, ctx);

    ca_sqr(C, D1, ctx);
    ca_sqr(t, D0, ctx);
    ca_mul(t, t, D0, ctx);
    ca_mul_ui(t, t, 4, ctx);
    ca_sub(C, C, t, ctx);
    ca_sqrt(C, C, ctx);

    /* C = (D1 + sqrt(D1^2 - 4*D0^3)) / 2  or
           (D1 - sqrt(D1^2 - 4*D0^3)) / 2,
            whichever is nonzero */
    success = 1;
    ca_add(t, D1, C, ctx);
    if (ca_check_is_zero(t, ctx) == T_FALSE)
    {
        ca_swap(C, t, ctx);
    }
    else if (ca_check_is_zero(t, ctx) != T_FALSE)
    {
        ca_sub(t, D1, C, ctx);
        if (ca_check_is_zero(t, ctx) == T_FALSE)
            ca_swap(C, t, ctx);
        else
            success = 0;
    }

    if (success)
    {
        ca_div_ui(C, C, 2, ctx);

        /* C = C^(1/3) */
        ca_set_ui(t, 1, ctx);
        ca_div_ui(t, t, 3, ctx);
        ca_pow(C, C, t, ctx);

        /* w1 = exp(2 pi i * (1 / 3)) */
        /* w2 = exp(2 pi i * (2 / 3)) = w1^2 */
        ca_omega(w1, ctx);
        ca_sqr(w2, w1, ctx);

        /* r1 = (b + C + D0 / C) / (-3*a) */
        /* r2 = (b + w1*C + w2 * D0 / C) / (-3*a) */
        /* r3 = (b + w2*C + w1 * D0 / C) / (-3*a) */

        ca_div(r1, D0, C, ctx);
        ca_mul(r2, r1, w2, ctx);
        ca_mul(r3, r1, w1, ctx);

        ca_add(r1, r1, C, ctx);
        ca_mul(t, w1, C, ctx);
        ca_add(r2, r2, t, ctx);
        ca_mul(t, w2, C, ctx);
        ca_add(r3, r3, t, ctx);

        ca_add(r1, r1, b, ctx);
        ca_add(r2, r2, b, ctx);
        ca_add(r3, r3, b, ctx);

        ca_mul_si(t, a, -3, ctx);
        ca_inv(t, t, ctx);

        ca_mul(r1, r1, t, ctx);
        ca_mul(r2, r2, t, ctx);
        ca_mul(r3, r3, t, ctx);
    }
    else
    {
        ca_unknown(r1, ctx);
        ca_unknown(r2, ctx);
        ca_unknown(r3, ctx);
    }

    ca_clear(D0, ctx);
    ca_clear(D1, ctx);
    ca_clear(C, ctx);
    ca_clear(w1, ctx);
    ca_clear(w2, ctx);
    ca_clear(t, ctx);

    return success;
}


int
_ca_poly_roots(ca_ptr roots, ca_srcptr poly, slong len, ca_ctx_t ctx)
{
    slong deg;
    truth_t leading_zero;

    if (len == 0)
        return 0;

    deg = len - 1;

    leading_zero = ca_check_is_zero(poly + deg, ctx);

    if (leading_zero != T_FALSE)
        return 0;

    while (deg > 1 && ca_check_is_zero(poly, ctx) == T_TRUE)
    {
        ca_zero(roots, ctx);
        roots++;
        poly++;
        len--;
        deg--;
    }

    if (deg == 0)
        return 1;

    if (deg == 1)
    {
        ca_div(roots, poly, poly + 1, ctx);
        ca_neg(roots, roots, ctx);
        return 1;
    }

    if (_ca_vec_is_fmpq_vec(poly, len, ctx))
    {
        fmpz_poly_t f;
        qqbar_ptr r;
        slong i;

        f->coeffs = _fmpz_vec_init(len);
        f->length = f->alloc = len;
        r = _qqbar_vec_init(len - 1);

        if (_ca_vec_fmpq_vec_is_fmpz_vec(poly, len, ctx))
        {
            for (i = 0; i < len; i++)
                fmpz_set(f->coeffs + i, CA_FMPQ_NUMREF(poly + i));
        }
        else
        {
            fmpz_t t;
            fmpz_init(t);
            fmpz_one(t);
            for (i = 0; i < len; i++)
                fmpz_lcm(t, t, CA_FMPQ_DENREF(poly + i));

            for (i = 0; i < len; i++)
            {
                fmpz_divexact(f->coeffs + i, t, CA_FMPQ_DENREF(poly + i));
                fmpz_mul(f->coeffs + i, f->coeffs + i, CA_FMPQ_NUMREF(poly + i));
            }

            fmpz_clear(t);
        }

        qqbar_roots_fmpz_poly(r, f, 0);

        for (i = 0; i < deg; i++)
            ca_set_qqbar(roots + i, r + i, ctx);

        _fmpz_vec_clear(f->coeffs, len);
        _qqbar_vec_clear(r, len - 1);

        return 1;
    }

    if (deg == 2)
    {
        _ca_poly_roots_quadratic(roots, roots + 1, poly + 2, poly + 1, poly + 0, ctx);
        return 1;
    }

    if (deg == 3)
    {
        return _ca_poly_roots_cubic(roots, roots + 1, roots + 2,
            poly + 3, poly + 2, poly + 1, poly + 0, ctx);
    }

    return 0;
}

int
ca_poly_roots(ca_vec_t roots, ulong * exp, const ca_poly_t poly, ca_ctx_t ctx)
{
    ca_poly_vec_t fac;
    ca_t c;
    slong i, j, num_roots, factor_deg;
    int success;
    ulong * fac_exp;

    if (poly->length == 0)
        return 0;

    ca_poly_vec_init(fac, 0, ctx);
    ca_init(c, ctx);
    fac_exp = flint_malloc(sizeof(ulong) * poly->length);
    success = ca_poly_factor_squarefree(c, fac, fac_exp, poly, ctx);

    if (success)
    {
        num_roots = 0;
        for (i = 0; i < fac->length; i++)
            num_roots += (fac->entries + i)->length - 1;

        ca_vec_set_length(roots, num_roots, ctx);

        num_roots = 0;
        for (i = 0; success && i < fac->length; i++)
        {
            factor_deg = (fac->entries + i)->length - 1;

            success = _ca_poly_roots(roots->entries + num_roots, (fac->entries + i)->coeffs, (fac->entries + i)->length, ctx);

            for (j = 0; j < factor_deg; j++)
                exp[num_roots + j] = fac_exp[i];

            num_roots += factor_deg;
        }
    }

    ca_poly_vec_clear(fac, ctx);
    ca_clear(c, ctx);
    flint_free(fac_exp);

    return success;
}
