/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2015 Vladimir Glazachev
   
******************************************************************************/

#include "aprcl.h"
#include <time.h>

int
unity_zp_is_inplace(ulong p, ulong exp)
{
    if (p == 2 && exp == 2) return 1;
    if (p == 2 && exp == 3) return 1;
    if (p == 2 && exp == 4) return 1;
    if (p == 3 && exp == 1) return 1;
    if (p == 3 && exp == 2) return 1;
    if (p == 5 && exp == 1) return 1;
    if (p == 7 && exp == 1) return 1;
    if (p == 11 && exp == 1) return 1;
    return 0;
}

void
unity_zp_pow_mont_fmpz(unity_zp f, const unity_zp g, const fmpz_t pow)
{
    fmpz_t r, ninv;
    unity_zp_mont fm, gm;

    if (unity_zp_is_inplace(f->p, f->exp) == 1)
    {
        unity_zp_pow_sliding_fmpz(f, g, pow);
        return;
    }

    fmpz_init(r);
    fmpz_init(ninv);
    unity_zp_mont_init(fm, f->p, f->exp, f->n, f->ninv);
    unity_zp_mont_init(gm, f->p, f->exp, f->n, f->ninv);
    unity_zp_mont_ninv(ninv, fm);
    fmpz_set(fm->ninv, ninv);
    fmpz_set(gm->ninv, ninv);

    fmpz_setbit(r, fm->r);

    fmpz_mul_2exp(fm->nr, fm->n, fm->r);
    fmpz_mul_2exp(gm->nr, gm->n, fm->r);

    unity_zp_to_mont(gm, g);
    _unity_zp_pow_mont_fmpz(fm, gm, pow, r);

    unity_zp_from_mont(f, fm);

    fmpz_clear(r);
    fmpz_clear(ninv);
    unity_zp_mont_clear(fm);
    unity_zp_mont_clear(gm);
}

/* Compute f = g^pow for f, g - unity_zp in Montgomery form */
void
_unity_zp_pow_mont_fmpz(unity_zp_mont f, const unity_zp_mont g,
        const fmpz_t pow, const fmpz_t r)
{
    ulong h, k, value;
    slong i, j;
    unity_zp_mont temp;
    unity_zp_mont *g_powers;

    fmpz_t * t;
    t = (fmpz_t*) flint_malloc(sizeof(fmpz_t) * (50));
    for (i = 0; i < 50; i++)
        fmpz_init(t[i]);


    unity_zp_mont_init(temp, f->p, f->exp, f->n, f->ninv);
    fmpz_set(temp->nr, f->nr);

    /* temp = g * g */
    unity_zp_mont_sqr(temp, g);

    /* selects optimal k value for n */
    k = _unity_zp_pow_2k_find_k(pow);

    /* 
        g_powers store odd powers of g up to 2^k - 1;
        g_powers[(i + 1) / 2] = g^i
    */
    g_powers = (unity_zp_mont*) flint_malloc(
            sizeof(unity_zp_mont) * (n_pow(2, k - 1) + 1));

    /* sets g_powers[0] = 1 */
    unity_zp_mont_init(g_powers[0], f->p, f->exp, f->n, f->ninv);
    fmpz_set(g_powers[0]->nr, f->nr);
    fmpz_poly_set_coeff_fmpz(g_powers[0]->poly, 0, r);

    /* sets g_powers[1] = g */
    unity_zp_mont_init(g_powers[1], f->p, f->exp, f->n, f->ninv);
    fmpz_set(g_powers[1]->nr, f->nr);
    unity_zp_mont_mul(g_powers[1], g_powers[0], g);

    /* sets g_powers[i] = g^2 * g_powers[i - 1] */
    for (i = 2; i <= n_pow(2, k - 1); i++)
    {
        unity_zp_mont_init(g_powers[i], f->p, f->exp, f->n, f->ninv);
        fmpz_set(g_powers[i]->nr, f->nr);
        unity_zp_mont_mul(g_powers[i], g_powers[i - 1], temp);
    }

    unity_zp_mont_set_zero(f);
    fmpz_poly_set_coeff_fmpz(f->poly, 0, r);
    i = fmpz_bits(pow) - 1;

    /* working with pow = (e_l, e_{l-1}, ... , e_0) in 2 base */
    while (i >= 0)
    {
        if (fmpz_tstbit(pow, i) == 0)
        {
            unity_zp_mont_sqr_inplace(temp, f, t);
            unity_zp_mont_swap(temp, f);
            i--;
        }
        else
        {
            /* 
                finds length of chain; chain is length of
                longest bitstring less then k ending on 1
            */
            j = FLINT_MAX(i - k + 1, 0);
            while (fmpz_tstbit(pow, j) == 0 && j <= i)
                j++;

            /* f = f^(2^(i - j + 1)) */
            for (h = 0; h < i - j + 1; h++)
            {
                unity_zp_mont_sqr_inplace(temp, f, t);
                unity_zp_mont_swap(temp, f);
            }

            /* 
                value = binary number (e_i, ... , e_j) in decimal base
            */
            value = 0;
            for (h = 0; h < i - j + 1; h++)
                value += fmpz_tstbit(pow, j + h) << h;

            /* f = f * g^value */
            unity_zp_mont_mul(temp, f, g_powers[(value + 1) / 2]);
            unity_zp_mont_swap(temp, f);

            /* increase i */
            i = j - 1;
        }
    }

    for (i = 0; i < 50; i++)
        fmpz_clear(t[i]);
    flint_free(t);

    for (i = 0; i <= n_pow(2, k - 1); i++)
        unity_zp_mont_clear(g_powers[i]);
    flint_free(g_powers);

    unity_zp_mont_clear(temp);
}

