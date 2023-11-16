/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_hypgeom.h"
#include "bernoulli.h"

static const unsigned int central_bin_tab[] = {
    1, 2, 6, 20, 70, 252, 924, 3432, 12870, 48620, 184756, 705432, 2704156,
    10400600, 40116600, 155117520, 601080390, 2333606220U,
};

void arb_gamma_stirling_coeff(arb_t b, ulong k, int digamma, slong prec);

/* See Richard P. Brent, "Asymptotic approximation of central binomial
coefficients with rigorous error bounds". https://arxiv.org/abs/1608.04834 */
static void
arb_hypgeom_central_bin_ui_asymp(arb_t res, ulong n, slong prec)
{
    arb_t s, t, u;
    fmpz_t n2;
    slong j, k, term_prec, wp;
    double term_mag, n2_mag;
    mag_t err, err2;

    arb_init(s);
    arb_init(t);
    arb_init(u);
    fmpz_init(n2);
    mag_init(err);
    mag_init(err2);

    wp = prec + 8;

    n2_mag = log(n) * 1.44269504088896;

    for (k = 1; k < prec; k++)
    {
        term_mag = bernoulli_bound_2exp_si(2 * k + 2) - (2 * k + 1) * n2_mag;
        term_mag -= (FLINT_BIT_COUNT((k + 1)*(2*k+1)) - 1);
        if (term_mag < -wp)
            break;
    }

    wp += 2 * FLINT_BIT_COUNT(k);

    BERNOULLI_ENSURE_CACHED(2*k)

    fmpz_set_ui(n2, n);
    fmpz_mul_ui(n2, n2, n);

    n2_mag *= 2;

    for (j = 0; j <= k - 1; j++)
    {
        term_mag = bernoulli_bound_2exp_si(2 * j + 2);
        term_mag -= j * n2_mag;
        term_prec = wp + term_mag;
        term_prec = FLINT_MIN(term_prec, wp);
        term_prec = FLINT_MAX(term_prec, 10);

        arb_gamma_stirling_coeff(t, j + 1, 0, term_prec);

        arb_mul_2exp_si(u, t, -2*j - 2);
        arb_sub(t, u, t, term_prec);
        arb_mul_2exp_si(t, t, 1);

        arb_addmul_fmpz(t, s, n2, wp);
        arb_swap(s, t);
    }

    arb_set_fmpz(t, n2);
    arb_pow_ui(t, t, k - 1, wp);
    arb_mul_ui(t, t, n, wp);
    arb_div(s, s, t, wp);

    /* error term: bernoulli(2k+2) / ((k+1)(2k+1)) / n^(2k+1) */
    mag_bernoulli_div_fac_ui(err, 2 * k + 2);
    mag_fac_ui(err2, 2 * k + 2);
    mag_mul(err, err, err2);
    mag_set_ui_lower(err2, n);
    mag_pow_ui_lower(err2, err2, 2 * k + 1);
    mag_mul_ui_lower(err2, err2, k + 1);
    mag_div(err, err, err2);
    arb_add_error_mag(s, err);

    arb_exp(s, s, wp);

    arb_const_pi(t, wp);
    arb_mul_ui(t, t, n, wp);
    arb_rsqrt(t, t, wp);

    arb_mul(res, s, t, prec);

    fmpz_set_ui(n2, n);
    fmpz_mul_2exp(n2, n2, 1);
    arb_mul_2exp_fmpz(res, res, n2);

    arb_clear(s);
    arb_clear(t);
    arb_clear(u);
    fmpz_clear(n2);
    mag_clear(err);
    mag_clear(err2);
}

void
arb_hypgeom_central_bin_ui(arb_t res, ulong n, slong prec)
{
    if (n <= 17)
    {
        arb_set_ui(res, central_bin_tab[n]);
        arb_set_round(res, res, prec);
    }
    else if (n < 6.0 * prec + 200.0)
    {
        fmpz_t t;
        fmpz_init(t);
        fmpz_bin_uiui(t, 2 * n, n);
        arb_set_round_fmpz(res, t, prec);
        fmpz_clear(t);
    }
    else
    {
        arb_hypgeom_central_bin_ui_asymp(res, n, prec);
    }
}

