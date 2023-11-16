/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_vec.h"
#include "arb.h"
#include "bernoulli.h"

void
bernoulli_rev_init(bernoulli_rev_t iter, ulong nmax)
{
    slong j;
    fmpz_t t;
    arb_t x;
    arf_t u;
    int round1, round2;
    slong wp;

    nmax -= (nmax % 2);
    iter->n = nmax;

    iter->alloc = 0;
    if (nmax < BERNOULLI_REV_MIN)
        return;

    iter->prec = wp = bernoulli_global_prec(nmax);

    iter->max_power = bernoulli_zeta_terms(nmax, iter->prec);
    iter->alloc = iter->max_power + 1;
    iter->powers = _fmpz_vec_init(iter->alloc);
    fmpz_init(iter->pow_error);
    arb_init(iter->prefactor);
    arb_init(iter->two_pi_squared);

    arb_init(x);
    fmpz_init(t);
    arf_init(u);

    /* precompute powers */
    for (j = 3; j <= iter->max_power; j += 2)
    {
        arb_ui_pow_ui(x, j, nmax, bernoulli_power_prec(j, nmax, wp));
        arb_inv(x, x, bernoulli_power_prec(j, nmax, wp));
        round1 = arf_get_fmpz_fixed_si(t, arb_midref(x), -wp);
        fmpz_set(iter->powers + j, t);

        /* error: the radius, plus two roundings */
        arf_set_mag(u, arb_radref(x));
        round2 = arf_get_fmpz_fixed_si(t, u, -wp);
        fmpz_add_ui(t, t, (round1 != 0) + (round2 != 0));
        if (fmpz_cmp(iter->pow_error, t) < 0)
            fmpz_set(iter->pow_error, t);
    }

    /* precompute (2pi)^2 and 2*(n!)/(2pi)^n */
    arb_fac_ui(iter->prefactor, nmax, wp);
    arb_mul_2exp_si(iter->prefactor, iter->prefactor, 1);

    arb_const_pi(x, wp);
    arb_mul_2exp_si(x, x, 1);
    arb_mul(iter->two_pi_squared, x, x, wp);

    arb_pow_ui(x, iter->two_pi_squared, nmax / 2, wp);
    arb_div(iter->prefactor, iter->prefactor, x, wp);

    fmpz_clear(t);
    arb_clear(x);
    arf_clear(u);
}

