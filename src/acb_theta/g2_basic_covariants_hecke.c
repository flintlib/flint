/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"
#include "profiler.h"

static void
hecke_coset(fmpz_mat_t m, slong k, slong p)
{
    slong a, b, c;
    slong i;

    if ((k < 0) || (k >= acb_theta_g2_hecke_nb(p)))
    {
        return;
    }

    fmpz_mat_zero(m);

    if (k < n_pow(p, 3))
    {
        /* Case 1 */
        a = k % p;
        b = (k / p) % p;
        c = (k / n_pow(p, 2)) % p;
        for (i = 0; i < 2; i++)
        {
            fmpz_one(fmpz_mat_entry(m, i, i));
        }
        for (i = 2; i < 4; i++)
        {
            fmpz_set_si(fmpz_mat_entry(m, i, i), p);
        }
        fmpz_set_si(fmpz_mat_entry(m, 0, 2), a);
        fmpz_set_si(fmpz_mat_entry(m, 0, 3), b);
        fmpz_set_si(fmpz_mat_entry(m, 1, 2), b);
        fmpz_set_si(fmpz_mat_entry(m, 1, 3), c);
    }
    else if (k < 1 + n_pow(p, 3))
    {
        /* Case 2 */
        fmpz_set_si(fmpz_mat_entry(m, 0, 0), p);
        fmpz_set_si(fmpz_mat_entry(m, 1, 1), p);
        fmpz_set_si(fmpz_mat_entry(m, 2, 2), 1);
        fmpz_set_si(fmpz_mat_entry(m, 3, 3), 1);
    }
    else if (k < 1 + n_pow(p, 3) + p)
    {
        /* Case 3 */
        a = k - n_pow(p, 3) - 1;
        fmpz_set_si(fmpz_mat_entry(m, 0, 0), 1);
        fmpz_set_si(fmpz_mat_entry(m, 0, 2), a);
        fmpz_set_si(fmpz_mat_entry(m, 1, 1), p);
        fmpz_set_si(fmpz_mat_entry(m, 2, 2), p);
        fmpz_set_si(fmpz_mat_entry(m, 3, 3), 1);
    }
    else
    {
        /* Case 4 */
        a = (k - 1 - n_pow(p, 3) - p) % p;
        b = ((k - 1 - n_pow(p, 3) - p)/p) % p;
        fmpz_set_si(fmpz_mat_entry(m, 0, 0), p);
        fmpz_set_si(fmpz_mat_entry(m, 1, 0), -a);
        fmpz_set_si(fmpz_mat_entry(m, 1, 1), 1);
        fmpz_set_si(fmpz_mat_entry(m, 1, 3), b);
        fmpz_set_si(fmpz_mat_entry(m, 2, 2), 1);
        fmpz_set_si(fmpz_mat_entry(m, 2, 3), a);
        fmpz_set_si(fmpz_mat_entry(m, 3, 3), p);
    }
}

static void
hecke_T1_coset(fmpz_mat_t m, slong k, slong p)
{
    slong a, b, c;
    slong i;

    if ((k < 0) || (k >= acb_theta_g2_hecke_nb(p * p)))
    {
        return;
    }

    fmpz_mat_zero(m);

    if (k == 0)
    {
        /* Case 1 */
        fmpz_set_si(fmpz_mat_entry(m, 0, 0), p);
        fmpz_set_si(fmpz_mat_entry(m, 1, 1), n_pow(p, 2));
        fmpz_set_si(fmpz_mat_entry(m, 2, 2), p);
        fmpz_set_si(fmpz_mat_entry(m, 3, 3), 1);
    }
    else if (k < 1 + (n_pow(p, 2)-1) )
    {
        /* Case 2 */
        if (k < 1 + (p-1))
        {
            /* a is zero, b too, c is anything nonzero */
            a = 0;
            b = 0;
            c = k;
        }
        else
        {
            /* a is nonzero, b is anything, c is b^2/a */
            /* k-p is between 0 and p(p-1)-1 */
            a = (k-p) % (p-1);
            a += 1;
            b = (k-p) % p;
            c = (b*b) % p;
            c *= n_invmod(a, p);
            c = c % p;
        }
        for (i = 0; i < 4; i++) fmpz_set_si(fmpz_mat_entry(m, i, i), p);
        fmpz_set_si(fmpz_mat_entry(m, 0, 2), a);
        fmpz_set_si(fmpz_mat_entry(m, 0, 3), b);
        fmpz_set_si(fmpz_mat_entry(m, 1, 2), b);
        fmpz_set_si(fmpz_mat_entry(m, 1, 3), c);
    }
    else if (k < n_pow(p, 2) + p)
    {
        /* Case 3 */
        a = k - n_pow(p, 2);
        fmpz_set_si(fmpz_mat_entry(m, 0, 0), n_pow(p, 2));
        fmpz_set_si(fmpz_mat_entry(m, 1, 0), -a*p);
        fmpz_set_si(fmpz_mat_entry(m, 1, 1), p);
        fmpz_set_si(fmpz_mat_entry(m, 2, 2), 1);
        fmpz_set_si(fmpz_mat_entry(m, 2, 3), a);
        fmpz_set_si(fmpz_mat_entry(m, 3, 3), p);
    }
    else if (k < n_pow(p, 2) + p + n_pow(p, 3))
    {
        /* Case 4 */
        k = k - n_pow(p, 2) - p;
        b = k % p;
        a = k / p;
        fmpz_set_si(fmpz_mat_entry(m, 0, 0), 1);
        fmpz_set_si(fmpz_mat_entry(m, 0, 2), a);
        fmpz_set_si(fmpz_mat_entry(m, 0, 3), -b);
        fmpz_set_si(fmpz_mat_entry(m, 1, 1), p);
        fmpz_set_si(fmpz_mat_entry(m, 1, 2), -p*b);
        fmpz_set_si(fmpz_mat_entry(m, 2, 2), n_pow(p, 2));
        fmpz_set_si(fmpz_mat_entry(m, 3, 3), p);
    }
    else
    {
        /* Case 5 */
        k = k - n_pow(p, 3) - n_pow(p, 2) - p;
        a = k%p;
        k = k/p;
        b = k%p;
        c = k/p;
        fmpz_set_si(fmpz_mat_entry(m, 0, 0), p);
        fmpz_set_si(fmpz_mat_entry(m, 0, 3), b*p);
        fmpz_set_si(fmpz_mat_entry(m, 1, 0), -a);
        fmpz_set_si(fmpz_mat_entry(m, 1, 1), 1);
        fmpz_set_si(fmpz_mat_entry(m, 1, 2), b);
        fmpz_set_si(fmpz_mat_entry(m, 1, 3), a*b+c);
        fmpz_set_si(fmpz_mat_entry(m, 2, 2), p);
        fmpz_set_si(fmpz_mat_entry(m, 2, 3), a*p);
        fmpz_set_si(fmpz_mat_entry(m, 3, 3), n_pow(p, 2));
    }
}

static void
hecke_T1_covariants(acb_ptr res, const acb_mat_t tau, slong p, slong prec)
{
    slong nb = ACB_THETA_G2_BASIC_NB;
    fmpz_mat_t mat;
    acb_mat_t w, c, cinv;
    acb_poly_t r;
    slong k;

    fmpz_mat_init(mat, 4, 4);
    acb_mat_init(w, 2, 2);
    acb_mat_init(c, 2, 2);
    acb_mat_init(cinv, 2, 2);
    acb_poly_init(r);

    for (k = 0; k < acb_theta_g2_hecke_nb(p * p); k++)
    {
        hecke_T1_coset(mat, k, p);
        acb_siegel_transform_cocycle_inv(w, c, cinv, mat, tau, prec);
        acb_theta_g2_fundamental_covariant(r, w, prec);
        acb_theta_g2_detk_symj(r, cinv, r, -2, 6, prec);
        acb_theta_g2_basic_covariants_lead(res + nb * k, r, prec);
    }

    fmpz_mat_clear(mat);
    acb_mat_clear(w);
    acb_mat_clear(c);
    acb_mat_clear(cinv);
    acb_poly_clear(r);
}

static void
hecke_covariants(acb_ptr res, const acb_mat_t tau, slong p, slong prec)
{
    slong nb = ACB_THETA_G2_BASIC_NB;
    fmpz_mat_t mat;
    acb_mat_t w, c, cinv;
    acb_poly_t r;
    slong k;

    fmpz_mat_init(mat, 4, 4);
    acb_mat_init(w, 2, 2);
    acb_mat_init(c, 2, 2);
    acb_mat_init(cinv, 2, 2);
    acb_poly_init(r);

    for (k = 0; k < acb_theta_g2_hecke_nb(p); k++)
    {
        hecke_coset(mat, k, p);
        acb_siegel_transform_cocycle_inv(w, c, cinv, mat, tau, prec);
        acb_theta_g2_fundamental_covariant(r, w, prec);
        acb_theta_g2_detk_symj(r, cinv, r, -2, 6, prec);
        acb_theta_g2_basic_covariants_lead(res + nb * k, r, prec);
    }

    fmpz_mat_clear(mat);
    acb_mat_clear(w);
    acb_mat_clear(c);
    acb_mat_clear(cinv);
    acb_poly_clear(r);
}

void
acb_theta_g2_basic_covariants_hecke(acb_ptr res, const acb_mat_t tau, slong q, slong prec)
{
    slong nb = ACB_THETA_G2_BASIC_NB;
    slong p;
    acb_poly_t r;
    int is_T1;

    if (n_is_prime(q))
    {
        p = q;
        is_T1 = 0;
    }
    else
    {
        p = n_sqrt(q);
        is_T1 = 1;
        if (p * p != q || !n_is_prime(p))
        {
            return;
        }
    }

    acb_poly_init(r);

    acb_theta_g2_fundamental_covariant(r, tau, prec);
    acb_theta_g2_basic_covariants_lead(res, r, prec);
    if (is_T1)
    {
        hecke_T1_covariants(res + nb, tau, p, prec);
    }
    else
    {
        hecke_covariants(res + nb, tau, p, prec);
    }

    acb_poly_clear(r);
}
