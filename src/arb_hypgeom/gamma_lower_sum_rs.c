/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "arb_hypgeom.h"

#ifdef __GNUC__
# define fabs __builtin_fabs
# define ldexp __builtin_ldexp
# define log __builtin_log
#else
# include <math.h>
#endif

static slong
exp_series_prec(slong k, double dz, double logdz, slong prec)
{
    double gain;

    if (prec <= 128)
        return prec;

    if (k <= dz + 5 || k <= 5)
        return prec;

    gain = (dz - k * logdz + k * (log(k) - 1.0)) * 1.4426950408889634;
    gain = FLINT_MAX(gain, 0);

    prec = prec - gain;
    prec = FLINT_MAX(prec, 32);
    return prec;
}

void
_arb_hypgeom_gamma_lower_sum_rs_1(arb_t res, ulong p, ulong q, const arb_t z, slong N, slong prec)
{
    slong m, j, k, jlen, jbot, wp;
    double dz, logdz;
    mp_limb_t c, chi, clo;
    arb_t s;
    arb_ptr zpow;
    mp_ptr cs;

    m = n_sqrt(N);
    m = FLINT_MAX(m, 2);
    k = N - 1;
    j = k % m;
    jlen = 0;
    jbot = j;
    c = 1;

    dz = arf_get_d(arb_midref(z), ARF_RND_UP);
    dz = fabs(dz);

    if (arf_cmpabs_2exp_si(arb_midref(z), prec) >= 0)
    {
        dz = ldexp(1.0, prec);
        logdz = ARF_EXP(arb_midref(z)) * log(2);
    }
    else if (arf_cmpabs_2exp_si(arb_midref(z), -32) >= 0)
    {
        logdz = log(dz);
    }
    else if (arf_cmpabs_2exp_si(arb_midref(z), -prec) <= 0)
    {
        logdz = -prec * log(2);
    }
    else
    {
        logdz = ARF_EXP(arb_midref(z)) * log(2);
    }

    arb_init(s);
    zpow = _arb_vec_init(m + 1);
    cs = flint_malloc(sizeof(mp_limb_t) * (m + 1));
    arb_mul_ui(zpow + m, z, q, prec);
    _arb_vec_set_powers(zpow, zpow + m, m + 1, prec);

    while (k >= 0)
    {
        if (k != 0)
        {
            /* Check if new coefficient will overflow limb */
            umul_ppmm(chi, clo, c, p + (k - 1) * q);

            if (chi != 0)
            {
                wp = exp_series_prec(k, dz, logdz, prec);

                /* Denominator will change, so evaluate current dot product */
                if (jlen != 0)
                {
                    arb_dot_ui(s, s, 0, zpow + jbot, 1, cs + jbot, 1, jlen, wp);
                    jlen = 0;
                }

                arb_div_ui(s, s, c, wp);
                c = 1;
            }
        }

        /* Update dot product */
        cs[j] = c;
        jlen++;
        jbot = j;

        if (k != 0)
        {
            c = c * (p + (k - 1) * q);

            /* Giant-step time. */
            if (j == 0)
            {
                wp = exp_series_prec(k, dz, logdz, prec);

                /* Evaluate current dot product */
                if (jlen != 0)
                {
                    arb_dot_ui(s, s, 0, zpow + jbot, 1, cs + jbot, 1, jlen, wp);
                    jlen = 0;
                }

                arb_mul(s, s, zpow + m, wp);
                j = m - 1;
            }
            else
            {
                j--;
            }
        }

        k--;
    }

    if (jlen != 0)
    {
        arb_dot_ui(s, s, 0, zpow + jbot, 1, cs + jbot, 1, jlen, prec);
        jlen = 0;
    }

    arb_div_ui(res, s, c, prec);

    _arb_vec_clear(zpow, m + 1);
    arb_clear(s);
    flint_free(cs);
}
