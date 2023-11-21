/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

void
acb_dirichlet_lerch_phi_direct(acb_t res, const acb_t z, const acb_t s, const acb_t a, slong prec)
{
    slong N, Nmax, wp, n;
    int a_real;
    acb_t negs, t, u, sum;
    mag_t C, S, zmag, tail_bound, tm, tol;

    if (!acb_is_finite(z) || !acb_is_finite(s) || !acb_is_finite(a))
    {
        acb_indeterminate(res);
        return;
    }

    if (acb_contains_int(a) && !arb_is_positive(acb_realref(a)))
    {
        if (!(acb_is_int(s) && arb_is_nonpositive(acb_realref(s))))
        {
            acb_indeterminate(res);
            return;
        }
    }

    acb_init(negs);
    acb_init(t);
    acb_init(u);
    acb_init(sum);

    acb_neg(negs, s);

    mag_init(C);
    mag_init(S);
    mag_init(zmag);
    mag_init(tail_bound);
    mag_init(tm);
    mag_init(tol);

    a_real = acb_is_real(a);
    wp = prec + 10;

    acb_get_mag(zmag, z);

    /* first term: 1/a^s */
    acb_pow(sum, a, negs, wp);

    acb_get_mag(tol, sum);
    mag_mul_2exp_si(tol, tol, -wp);

    if (a_real)
    {
        /* Tail bound |z|^N / |(a+N)^s| * sum C^k, C = |z| * exp(max(0, -re(s)) / (a+N))  */
        arb_nonnegative_part(acb_realref(t), acb_realref(negs));
        arb_get_mag(S, acb_realref(t));
    }
    else
    {
        /* Tail bound |z|^N / |(a+N)^s| * sum C^k, C = |z| * exp(|s / (a+N)|)  */
        acb_get_mag(S, s);
    }

    Nmax = 100 * prec + 0.1 * prec * n_sqrt(prec);
    Nmax = FLINT_MAX(Nmax, 1);
    Nmax = FLINT_MIN(Nmax, WORD_MAX / 2);

    mag_inf(tail_bound);

    for (N = 1; N <= Nmax; N = FLINT_MAX(N+4, N*1.1))
    {
        acb_add_ui(t, a, N, 53);

        if (arb_is_positive(acb_realref(t)))
        {
            acb_get_mag_lower(C, t);
            mag_div(C, S, C);
            mag_exp(C, C);
            mag_mul(C, C, zmag);
            mag_geom_series(C, C, 0);

            if (mag_is_finite(C))
            {
                mag_pow_ui(tail_bound, zmag, N);
                mag_mul(tail_bound, tail_bound, C);
                acb_pow(t, t, negs, 53);
                acb_get_mag(C, t);
                mag_mul(tail_bound, tail_bound, C);

                if (mag_cmp(tail_bound, tol) <= 0)
                    break;
            }
            else
            {
                mag_inf(tail_bound);
            }
        }
    }

    if (mag_is_finite(tail_bound))
    {
        acb_one(t);

        for (n = 1; n < N; n++)
        {
            if (n % 8 == 0 && !acb_is_real(z))
                acb_pow_ui(t, z, n, wp);
            else
                acb_mul(t, t, z, wp);

            acb_add_ui(u, a, n, wp);
            acb_pow(u, u, negs, wp);
            acb_mul(u, t, u, wp);
            acb_add(sum, sum, u, wp);
        }

        if (acb_is_real(z) && acb_is_real(s) && acb_is_real(a))
            arb_add_error_mag(acb_realref(sum), tail_bound);
        else
            acb_add_error_mag(sum, tail_bound);

        acb_set_round(res, sum, prec);
    }
    else
    {
        acb_indeterminate(res);
    }

    mag_clear(C);
    mag_clear(S);
    mag_clear(zmag);
    mag_clear(tail_bound);
    mag_clear(tm);
    mag_clear(tol);

    acb_clear(negs);
    acb_clear(t);
    acb_clear(u);
    acb_clear(sum);
}
