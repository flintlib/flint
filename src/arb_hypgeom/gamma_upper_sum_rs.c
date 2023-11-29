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
# define log __builtin_log
#else
# include <math.h>
#endif

static slong
asymp_prec(slong k, double logdz, slong prec)
{
    double gain;

    if (prec <= 128)
        return prec;

    if (k <= 5)
        return prec;

    gain = (k * logdz - (k * (log(k) - 1.0))) * 1.4426950408889634 - 4;
    gain = FLINT_MAX(gain, 0);

    prec = prec - gain;
    prec = FLINT_MAX(prec, 32);
    return prec;
}

void
_arb_hypgeom_gamma_upper_sum_rs_1(arb_t res, ulong p, ulong q, const arb_t z, slong N, slong prec)
{
    slong m, i, j, k, jlen, jbot, jtop, wp;
    double dz, logdz;
    mp_limb_t c, chi, clo;
    arb_t s, t;
    arb_ptr zpow;
    mp_ptr cs;

    m = n_sqrt(N);
    m = FLINT_MAX(m, 2);
    k = N - 1;
    j = k % m;
    jlen = 0;
    jbot = j;

    if (arf_cmpabs_2exp_si(arb_midref(z), prec) >= 0)
    {
        logdz = ARF_EXP(arb_midref(z)) * log(2);
    }
    else if (arf_cmpabs_2exp_si(arb_midref(z), -32) >= 0)
    {
        dz = arf_get_d(arb_midref(z), ARF_RND_UP);
        dz = fabs(dz);
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
    arb_init(t);
    zpow = _arb_vec_init(m + 1);
    cs = flint_malloc(sizeof(mp_limb_t) * (m + 1));
    arb_mul_ui(zpow + m, z, q, prec);
    arb_inv(zpow + m, zpow + m, prec);
    _arb_vec_set_powers(zpow, zpow + m, m + 1, prec);

    while (k >= 0)
    {
        /* Find run of coefficients whose product fits in a limb */
        jlen = 1;
        jtop = jbot = k;

        if (jtop > 0)
        {
            c = p + q * (jtop - 1);
            while (jlen <= j)
            {
                if (jbot >= 2)
                {
                    umul_ppmm(chi, clo, c, p + q * (jbot - 2));

                    if (chi != 0)
                        break;

                    c = clo;
                }

                jbot--;
                jlen++;
            }
        }

        if (jbot != jtop - jlen + 1)
            flint_throw(FLINT_ERROR, "(%s)\n", __func__);

        /* Factors between jbot and jtop inclusive */
        if (jbot == 0)
            cs[0] = 1;
        else
            cs[0] = p + q * (jbot - 1);

        for (i = 1; i < jlen; i++)
            cs[i] = cs[i - 1] * (p + q * (jbot + i - 1));

        wp = asymp_prec(k - jlen, logdz, prec);

        /* todo: special case jlen == 1 */
        arb_add(t, s, zpow + j, wp);
        arb_swap(zpow + j, t);
        arb_dot_ui(s, NULL, 0, zpow + j - jlen + 1, 1, cs, 1, jlen, wp);
        arb_swap(zpow + j, t);

        k -= jlen;
        j -= (jlen - 1);

        if (j == 0 && k >= 1)
        {
            arb_mul(s, s, zpow + m, wp);
            j = m - 1;
        }
        else
        {
            j--;
        }
    }

    arb_swap(res, s);

    _arb_vec_clear(zpow, m + 1);
    arb_clear(s);
    arb_clear(t);
    flint_free(cs);
}
