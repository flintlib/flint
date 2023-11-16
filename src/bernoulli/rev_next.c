/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "bernoulli.h"

static inline void
mag_ui_div(mag_t z, ulong c, const mag_t x)
{
    mag_t t;
    mag_init(t);
    mag_set_ui(t, c);
    mag_div(z, t, x);
    mag_clear(t);
}

static const short small_denom[] = {
    1, 6, 30, 42, 30, 66, 2730, 6, 510, 798, 330, 138, 2730, 6, 870, 14322
};

static const int small_numer[] = {
    1, 1, -1, 1, -1, 5, -691, 7, -3617, 43867, -174611, 854513, -236364091, 8553103
};

static void
_fmpq_bernoulli_small(fmpz_t p, fmpz_t q, ulong n)
{
    if (n == 1)
    {
        fmpz_set_si(p, -1);
        fmpz_set_ui(q, 2);
    }
    else if (n % 2 == 1)
    {
        fmpz_zero(p);
        fmpz_one(q);
    }
    else
    {
        if (n == 28)
            fmpz_set_d(p, -23749461029.0);
        else if (n == 30)
            fmpz_set_d(p, 8615841276005.0);
        else
            fmpz_set_si(p, small_numer[n / 2]);

        fmpz_set_si(q, small_denom[n / 2]);
    }
}

void
bernoulli_rev_next(fmpz_t numer, fmpz_t denom, bernoulli_rev_t iter)
{
    ulong n;
    slong j, wp;
    fmpz_t sum;
    mag_t err;
    arb_t z, h;

    n = iter->n;
    wp = iter->prec;

    if (n < BERNOULLI_REV_MIN)
    {
        _fmpq_bernoulli_small(numer, denom, n);
        if (n != 0)
            iter->n -= 2;
        return;
    }

    fmpz_init(sum);
    mag_init(err);
    arb_init(z);
    arb_init(h);

    /* add all odd powers */
    fmpz_zero(sum);
    for (j = iter->max_power; j >= 3; j -= 2)
        fmpz_add(sum, sum, iter->powers + j);
    arb_set_fmpz(z, sum);

    /* bound numerical error from the powers */
    fmpz_mul_ui(sum, iter->pow_error, iter->max_power / 2);
    mag_set_fmpz(err, sum);
    mag_add(arb_radref(z), arb_radref(z), err);
    arb_mul_2exp_si(z, z, -wp);

    arb_add_ui(z, z, 1, wp);

    /* add truncation error: sum_{k > N} 1/k^n <= 1/N^(i-1) */
    mag_set_ui_lower(err, iter->max_power);
    mag_pow_ui_lower(err, err, n - 1);
    mag_ui_div(err, 1, err);
    mag_add(arb_radref(z), arb_radref(z), err);

    /* convert zeta to Bernoulli number */
    arb_div_2expm1_ui(h, z, n, wp);
    arb_add(z, z, h, wp);
    arb_mul(z, z, iter->prefactor, wp);
    arith_bernoulli_number_denom(denom, n);
    arb_mul_fmpz(z, z, denom, wp);
    if (n % 4 == 0)
        arb_neg(z, z);

    /* flint_printf("%wd: ", n); arb_printd(z, 5); flint_printf("\n"); */

    if (!arb_get_unique_fmpz(numer, z))
    {
        flint_printf("warning: insufficient precision for B_%wd\n", n);
        _bernoulli_fmpq_ui(numer, denom, n);
    }

    /* update prefactor */
    if (n > 0)
    {
        arb_mul(iter->prefactor, iter->prefactor, iter->two_pi_squared, wp);
        arb_div_ui(iter->prefactor, iter->prefactor, n, wp);
        arb_div_ui(iter->prefactor, iter->prefactor, n - 1, wp);
    }

    /* update powers */
    for (j = 3; j <= iter->max_power; j += 2)
        fmpz_mul2_uiui(iter->powers + j, iter->powers + j, j, j);

    /* bound error after update */
    fmpz_mul2_uiui(iter->pow_error, iter->pow_error,
        iter->max_power, iter->max_power);

    /* readjust precision */
    if (n % 64 == 0 && n > BERNOULLI_REV_MIN)
    {
        slong new_prec, new_max;

        new_prec = bernoulli_global_prec(n);
        new_max = bernoulli_zeta_terms(n, new_prec);

        if (new_prec < iter->prec && new_max <= iter->max_power)
        {
            /* change precision of the powers */
            for (j = 3; j <= new_max; j += 2)
                fmpz_tdiv_q_2exp(iter->powers + j, iter->powers + j,
                    iter->prec - new_prec);

            /* the error also changes precision */
            fmpz_cdiv_q_2exp(iter->pow_error, iter->pow_error,
                iter->prec - new_prec);
            /* contribution of rounding error when changing the precision
               of the powers */
            fmpz_add_ui(iter->pow_error, iter->pow_error, 1);

            /* speed improvement (could be skipped with better multiplication) */
            arb_set_round(iter->two_pi_squared, iter->two_pi_squared, new_prec);

            iter->max_power = new_max;
            iter->prec = new_prec;
        }
    }

    iter->n -= 2;

    fmpz_clear(sum);
    mag_clear(err);
    arb_clear(z);
    arb_clear(h);
}

