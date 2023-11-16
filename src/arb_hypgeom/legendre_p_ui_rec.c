/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_hypgeom.h"
#include "gmpcompat.h"

void
arb_hypgeom_legendre_p_ui_rec(arb_t res, arb_t res_prime, ulong n, const arb_t x, slong prec)
{
    slong wp;
    ulong k, den;
    mp_limb_t denlo, denhi;
    mpz_t p0, p1, xx, tt;
    fmpz_t fxx;
    int error;
    arb_t t, u, x2sub1;
    mag_t err1, err2, xrad;

    if (n > (UWORD(1) << (FLINT_BITS / 2 - 1)))
    {
        if (res != NULL) arb_indeterminate(res);
        if (res_prime != NULL) arb_indeterminate(res);
        return;
    }

    if (n == 0)
    {
        if (res != NULL) arb_one(res);
        if (res_prime != NULL) arb_zero(res_prime);
        return;
    }

    mag_init(xrad);

    arb_get_mag(xrad, x);
    /* error analysis assumes |x| < 1 */
    if (mag_cmp_2exp_si(xrad, 0) >= 0)
    {
        arb_hypgeom_legendre_p_ui_one(res, res_prime, n, x, n + 1, prec);
        mag_clear(xrad);
        return;
    }

    mpz_init(p0);
    mpz_init(p1);
    mpz_init(xx);
    mpz_init(tt);
    fmpz_init(fxx);
    arb_init(t);
    arb_init(u);
    arb_init(x2sub1);
    mag_init(err1);
    mag_init(err2);

    wp = -arf_abs_bound_lt_2exp_si(arb_midref(x));
    wp = FLINT_MAX(wp, 0);
    wp = FLINT_MIN(wp, prec);
    wp += prec + 2 * FLINT_BIT_COUNT(n + 2); /* (n+2)^2 >= 0.75(n+1)(n+2)+1 */

    arb_mul(x2sub1, x, x, ARF_PREC_EXACT);
    arb_neg(x2sub1, x2sub1);
    arb_add_ui(x2sub1, x2sub1, 1, wp);

    error = arf_get_fmpz_fixed_si(fxx, arb_midref(x), -wp);
    fmpz_get_mpz(xx, fxx);

    mag_set(xrad, arb_radref(x));
    if (error)
        mag_add_ui_2exp_si(xrad, xrad, 1, -wp);

    mpz_set_ui(p0, 1);
    mpz_mul_2exp(p0, p0, wp);
    mpz_set(p1, xx);

    den = 1;
    for (k = 1; k < n; k++)
    {
        mpz_mul(tt, p1, xx);
        mpz_tdiv_q_2exp(tt, tt, wp);
        flint_mpz_mul_ui(p0, p0, k*k);
        mpz_neg(p0, p0);
        flint_mpz_addmul_ui(p0, tt, 2 * k + 1);
        mpz_swap(p0, p1);
        umul_ppmm(denhi, denlo, den, k + 1);
        if (denhi != 0)
        {
            flint_mpz_tdiv_q_ui(p0, p0, den);
            flint_mpz_tdiv_q_ui(p1, p1, den);
            den = k + 1;
        }
        else
        {
            den = denlo;
        }
    }
    flint_mpz_tdiv_q_ui(p0, p0, den/n);
    flint_mpz_tdiv_q_ui(p1, p1, den);

    if (!mag_is_zero(xrad))
    {
        arb_hypgeom_legendre_p_ui_deriv_bound(err1, err2, n, x, x2sub1);
        mag_mul(err1, err1, xrad);
        mag_mul(err2, err2, xrad);
    }

    arf_set_mpz(arb_midref(t), p1);
    arf_mul_2exp_si(arb_midref(t), arb_midref(t), -wp);
    mag_set_ui_2exp_si(arb_radref(t), (n + 1) * (n + 2), -wp);
    mag_add(arb_radref(t), arb_radref(t), err1);

    if (res_prime != NULL)
    {
        /* P' = n (P[n-1] - x P) / (1 - x^2) */
        arf_set_mpz(arb_midref(u), p0);
        arf_mul_2exp_si(arb_midref(u), arb_midref(u), -wp);
        mag_set_ui_2exp_si(arb_radref(u), n * (n + 1), -wp);

        arb_submul(u, t, x, wp);
        arb_div(res_prime, u, x2sub1, wp);
        arb_mul_ui(res_prime, res_prime, n, prec);

        mag_add(arb_radref(res_prime), arb_radref(res_prime), err2);
    }

    if (res != NULL)  /* done last since x may be aliased with res */
    {
        arb_set_round(res, t, prec);
    }

    mpz_clear(p0);
    mpz_clear(p1);
    mpz_clear(xx);
    mpz_clear(tt);
    fmpz_clear(fxx);
    arb_clear(t);
    arb_clear(u);
    arb_clear(x2sub1);
    mag_clear(err1);
    mag_clear(err2);
    mag_clear(xrad);
}

