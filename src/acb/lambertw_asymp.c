/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_vec.h"
#include "acb.h"

void
acb_lambertw_asymp(acb_t res, const acb_t z, const fmpz_t k, slong L, slong M, slong prec)
{
    acb_t L1, L2, sigma, tau, s, c, u;
    slong l, m;
    fmpz_t t;
    fmpz * sc;

    /* For k = 0, the asymptotic expansion is not valid near 0. */
    /* (It is sufficient to look at the midpoint as a test here.) */
    if (fmpz_is_zero(k) && arf_cmpabs_2exp_si(arb_midref(acb_realref(z)), 0) < 0
                        && arf_cmpabs_2exp_si(arb_midref(acb_imagref(z)), 0) < 0)
    {
        acb_indeterminate(res);
        return;
    }

    acb_init(L1);
    acb_init(L2);
    acb_init(sigma);
    acb_init(tau);
    acb_init(s);
    acb_init(c);
    acb_init(u);
    fmpz_init(t);

    acb_const_pi(L2, prec);
    acb_mul_2exp_si(L2, L2, 1);
    acb_mul_fmpz(L2, L2, k, prec);
    acb_mul_onei(L2, L2);
    acb_log(L1, z, prec);
    acb_add(L1, L1, L2, prec);
    acb_log(L2, L1, prec);

    acb_inv(sigma, L1, prec);
    acb_mul(tau, L2, sigma, prec);

    acb_zero(s);

    /* Stirling numbers */
    sc = _fmpz_vec_init(L);

    acb_one(u);

    for (m = 1; m < M; m++)
    {
        if (m == 1)
        {
            for (l = 0; l < L; l++)
                fmpz_one(sc + l);
        }
        else
        {
            for (l = 0; l < L; l++)
            {
                fmpz_mul_ui(sc + l, sc + l, m + l - 1);
                if (l > 0)
                    fmpz_add(sc + l, sc + l, sc + l - 1);
            }
        }

        acb_zero(c);

        /* todo: precompute powers instead of horner */
        for (l = L - 1; l >= 0; l--)
        {
            acb_mul(c, c, sigma, prec);
            if (l % 2)
                acb_sub_fmpz(c, c, sc + l, prec);
            else
                acb_add_fmpz(c, c, sc + l, prec);
        }

        acb_mul(u, u, tau, prec);
        acb_div_ui(u, u, m, prec);
        acb_addmul(s, c, u, prec);
    }

    _fmpz_vec_clear(sc, L);

    acb_sub(s, s, L2, prec);
    acb_add(s, s, L1, prec);

    {
        mag_t m4s, m4t, one, q, r;

        mag_init(m4s);
        mag_init(m4t);
        mag_init(one);
        mag_init(q);
        mag_init(r);

        acb_get_mag(m4s, sigma);
        mag_mul_2exp_si(m4s, m4s, 2);
        acb_get_mag(m4t, tau);
        mag_mul_2exp_si(m4t, m4t, 2);

        mag_one(one);

        mag_sub_lower(q, one, m4s);
        mag_sub_lower(r, one, m4t);
        mag_mul(q, q, r);

        mag_pow_ui(r, m4s, L);
        mag_mul(r, r, m4t);
        mag_pow_ui(m4t, m4t, M);
        mag_add(r, r, m4t);

        mag_div(q, r, q);

        acb_add_error_mag(s, q);

        mag_clear(m4s);
        mag_clear(m4t);
        mag_clear(one);
        mag_clear(q);
        mag_clear(r);
    }

    acb_set(res, s);

    acb_clear(sigma);
    acb_clear(tau);
    acb_clear(s);
    acb_clear(c);
    acb_clear(L1);
    acb_clear(L2);
    acb_clear(u);
    fmpz_clear(t);
}

