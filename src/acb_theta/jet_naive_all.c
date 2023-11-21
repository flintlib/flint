/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"
#include "acb_modular.h"
#include "acb_theta.h"

/* Use a big ellipsoid to avoid complicated formulas for derivatives; this
   introduces powers of i in worker */

static void
worker(acb_ptr dth, acb_srcptr v1, acb_srcptr v2, const slong * precs, slong len,
    const acb_t cofactor, const slong * coords, slong ord, slong g, slong prec, slong fullprec)
{
    slong n = 1 << g;
    slong nb = acb_theta_jet_nb(ord, g);
    slong * tups;
    slong a0, a1;
    slong * dots;
    acb_ptr v3, aux;
    acb_t x, y;
    fmpz_t num, t;
    slong j, i;
    ulong b;

    tups = flint_malloc(g * nb * sizeof(slong));
    dots = flint_malloc(n * sizeof(slong));
    v3 = _acb_vec_init(len);
    aux = _acb_vec_init(nb * n * n);
    acb_init(x);
    acb_init(y);
    fmpz_init(num);
    fmpz_init(t);

    /* Precompute a0, a1, dots */
    a0 = acb_theta_char_get_a(coords, g);
    a1 = a0 ^ (1 << (g - 1));
    for (b = 0; b < n; b++)
    {
        dots[b] = acb_theta_char_dot_slong(b, coords, g);
    }

    /* Compute products in v3 */
    for (i = 0; i < len; i++)
    {
        acb_mul(&v3[i], &v1[i], &v2[i], precs[i]);
    }

    acb_theta_jet_tuples(tups, ord, g);
    for (j = 0; j < nb; j++)
    {
        fmpz_one(num);
        for (i = 1; i < g; i++)
        {
            fmpz_set_si(t, coords[i]);
            fmpz_pow_ui(t, t, tups[j * g + i]);
            fmpz_mul(num, num, t);
        }

        /* Loop over lattice points */
        for (i = 0; i < len; i++)
        {
            fmpz_set_si(t, coords[0] + i);
            fmpz_pow_ui(t, t, tups[j * g]);
            acb_mul_fmpz(x, &v3[i], t, precs[i]);
            /* Loop over b, adding coefficients in both a0b and a1b */
            for (b = 0; b < n; b++)
            {
                acb_mul_i_pow_si(y, x, (dots[b] + i * (b >> (g - 1))) % 4);
                if (i % 2 == 0)
                {
                    acb_add(&aux[(n * a0 + b) * nb + j],
                        &aux[(n * a0 + b) * nb + j], y, prec);
                }
                else
                {
                    acb_add(&aux[(n * a1 + b) * nb + j],
                        &aux[(n * a1 + b) * nb + j], y, prec);
                }
            }
        }

        /* Multiply by cofactor * num */
        acb_mul_fmpz(x, cofactor, num, prec);
        for (b = 0; b < n; b++)
        {
            acb_mul(&aux[(n * a0 + b) * nb + j], &aux[(n * a0 + b) * nb + j], x, prec);
            acb_mul(&aux[(n * a1 + b) * nb + j], &aux[(n * a1 + b) * nb + j], x, prec);
        }
    }

    _acb_vec_add(dth, dth, aux, nb * n * n, fullprec);

    flint_free(tups);
    flint_free(dots);
    _acb_vec_clear(v3, len);
    _acb_vec_clear(aux, nb * n * n);
    acb_clear(x);
    acb_clear(y);
    fmpz_clear(num);
    fmpz_clear(t);
}

static void
acb_theta_jet_naive_all_gen(acb_ptr dth, acb_srcptr z, const acb_mat_t tau,
    slong ord, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n2 = 1 << (2 * g);
    slong nb = acb_theta_jet_nb(ord, g);
    slong * tups;
    acb_theta_eld_t E;
    arb_mat_t C;
    arf_t R2, eps;
    acb_ptr aux, new_z;
    acb_mat_t new_tau;
    arb_ptr v, a;
    acb_t c;
    arb_t u;
    fmpz_t m, t;
    slong k, j;
    int b;

    tups = flint_malloc(g * nb * sizeof(slong));
    acb_theta_eld_init(E, g, g);
    arb_mat_init(C, g, g);
    arf_init(R2);
    arf_init(eps);
    aux = _acb_vec_init(n2 * nb);
    new_z = _acb_vec_init(g);
    acb_mat_init(new_tau, g, g);
    v = _arb_vec_init(g);
    a = _arb_vec_init(g);
    acb_init(c);
    arb_init(u);
    fmpz_init(m);
    fmpz_init(t);

    _acb_vec_scalar_mul_2exp_si(new_z, z, g, -1);
    acb_mat_scalar_mul_2exp_si(new_tau, tau, -2);
    acb_siegel_cho(C, new_tau, prec);

    acb_theta_naive_reduce(v, new_z, a, c, u, new_z, 1, new_tau, prec);
    acb_theta_jet_naive_radius(R2, eps, C, v, ord, prec);
    b = acb_theta_eld_set(E, C, R2, v);

    if (b)
    {
        acb_theta_naive_worker(dth, nb * n2, new_z, 1, new_tau, E, ord, prec, worker);
        arb_mul_arf(u, u, eps, prec);
        for (k = 0; k < nb * n2; k++)
        {
            acb_mul(&dth[k], &dth[k], c, prec);
            acb_add_error_arb(&dth[k], u);
        }

        acb_theta_jet_tuples(tups, ord, g);
        for (k = 0; k < nb; k++)
        {
            acb_const_pi(c, prec); /* not 2 pi because of rescaling */
            acb_mul_onei(c, c);
            acb_pow_ui(c, c, acb_theta_jet_total_order(tups + k * g, g), prec);
            fmpz_one(m);
            for (j = 0; j < g; j++)
            {
                fmpz_fac_ui(t, tups[k * g + j]);
                fmpz_mul(m, m, t);
            }
            acb_div_fmpz(c, c, m, prec);
            for (j = 0; j < n2; j++)
            {
                acb_mul(&dth[j * nb + k], &dth[j * nb + k], c, prec);
            }
        }

        _arb_vec_neg(a, a, g);
        acb_theta_jet_exp_pi_i(aux, a, ord, g, prec);
        for (k = 0; k < n2; k++)
        {
            acb_theta_jet_mul(dth + k * nb, dth + k * nb, aux, ord, g, prec);
            arb_zero(u);
            for (j = 0; j < g; j++)
            {
                if ((k >> (g - 1 - j)) % 2 == 1)
                {
                    arb_add(u, u, &a[j], prec);
                }
            }
            acb_onei(c);
            acb_pow_arb(c, c, u, prec);
            _acb_vec_scalar_mul(dth + k * nb, dth + k * nb, nb, c, prec);
        }
    }
    else
    {
        for (k = 0; k < nb * n2; k++)
        {
            acb_indeterminate(&dth[k]);
        }
    }

    flint_free(tups);
    acb_theta_eld_clear(E);
    arb_mat_clear(C);
    arf_clear(R2);
    arf_clear(eps);
    _acb_vec_clear(aux, n2 * nb);
    _acb_vec_clear(new_z, g);
    acb_mat_clear(new_tau);
    _arb_vec_clear(v, g);
    _arb_vec_clear(a, g);
    acb_clear(c);
    arb_clear(u);
    fmpz_clear(m);
    fmpz_clear(t);
}

void
acb_theta_jet_naive_all(acb_ptr dth, acb_srcptr z, const acb_mat_t tau,
    slong ord, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong nb = acb_theta_jet_nb(ord, g);

    if (g == 1)
    {
        acb_modular_theta_jet(dth + 3 * nb, dth + 2 * nb, dth, dth + nb,
            z, acb_mat_entry(tau, 0, 0), nb, prec);
        _acb_vec_neg(dth + 3 * nb, dth + 3 * nb, nb);
    }
    else
    {
        acb_theta_jet_naive_all_gen(dth, z, tau, ord, prec);
    }
}
