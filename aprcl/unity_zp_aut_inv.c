/*
    Copyright (C) 2015 Vladimir Glazachev

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "aprcl.h"

/*
    Computes f such that \sigma_x(f) = g.
*/
void
unity_zp_aut_inv(unity_zp f, const unity_zp g, ulong x)
{
    ulong i, j, p_pow1, p_pow2, m, p_pow_preinv;
    fmpz_t f_coeff, g_coeff;
    
    fmpz_init(f_coeff);
    fmpz_init(g_coeff);
    p_pow1 = n_pow(f->p, f->exp - 1);   /* p_pow1 = p^{k - 1}       */
    p_pow2 = p_pow1 * f->p;             /* p_pow2 = p^k             */
    m = (f->p - 1) * p_pow1;            /* m = (p - 1) * p^{k - 1}  */
    p_pow_preinv = n_preinvert_limb(p_pow2);
    unity_zp_set_zero(f);

    /* for i = 0, 1,..., m - 1 set f[i] = g[xi mod p^k] */
    for (i = 0; i < m; i++)
    {
        /* set g_ind = x * i mod p^k */
        ulong g_ind = n_mulmod2_preinv(x, i, p_pow2, p_pow_preinv);

        /* set g_coeff to g[g_ind] */
        fmpz_mod_poly_get_coeff_fmpz(g_coeff, g->poly, g_ind, g->ctx);

        /* set f[i] = g[x * i mod p^k] */
        unity_zp_coeff_set_fmpz(f, i, g_coeff);
    }

    /*
        for i = m, m + 1,..., p^k - 1
        for j = 1, 2,..., p - 1 
        set f[i - j * p^{k - 1}] = 
        (f[i - j * p^{k - 1}] - g[x * i mod p^k]) mod n
    */
    for (i = m; i < p_pow2; i++)
    {
        /* set g_ind = x * i mod p^k */
        ulong g_ind = n_mulmod2_preinv(x, i, p_pow2, p_pow_preinv);

        for (j = 1; j < f->p; j++)
        {
            /* set f_ind = i - j * p^{k - 1} */
            ulong f_ind = i - j * p_pow1;

            /* set g_coeff = g[x * i mod p^k] */
            fmpz_mod_poly_get_coeff_fmpz(g_coeff, g->poly, g_ind, g->ctx);

            /* set f_coeff = f[i - j * p^{k - 1}] */
            fmpz_mod_poly_get_coeff_fmpz(f_coeff, f->poly, f_ind, f->ctx);

            /* set f_coeff = f[i - j * p^{k - 1}] - g[x * i mod p^k] */
            fmpz_sub(f_coeff, f_coeff, g_coeff);

            /* set f[i - j * p^{k - 1}] = f_coeff */
            unity_zp_coeff_set_fmpz(f, f_ind, f_coeff);
        }
    }

    fmpz_clear(f_coeff);
    fmpz_clear(g_coeff);
}

