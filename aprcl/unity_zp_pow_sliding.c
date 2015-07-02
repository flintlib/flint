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

void
unity_zp_pow_sliding_fmpz(unity_zp f, const unity_zp g, const fmpz_t pow)
{
    ulong h, k, value, pow2k, chain;
    ulong *digits;
    slong i, y, j;
    fmpz_t power;
    unity_zp g_sqr, temp;
    unity_zp *g_powers;

    fmpz_init_set(power, pow);
    unity_zp_init(g_sqr, f->p, f->exp, f->n);
    unity_zp_init(temp, f->p, f->exp, f->n);

    /* g_sqr = g * g */
    unity_zp_sqr(g_sqr, g);

    /* selects optimal k value for n */
    k = _unity_zp_pow_2k_find_k(pow);
    /* computes 2^k */
    pow2k = 2;

    /*
        set digits[i] = i-th digit of pow in binary representation;
        least signifant digit have number 0
    */
    digits = (ulong*) flint_malloc(sizeof(ulong) * fmpz_bits(pow));
    for (i = 0; i < fmpz_bits(pow); i++)
    {
        digits[i] = fmpz_tdiv_ui(power, 2);
        fmpz_tdiv_q_ui(power, power, 2);
    }

    /* 
        g_powers store odd powers of g up to 2^k - 1;
        g_powers[(i + 1) / 2] = g^i
    */
    g_powers = (unity_zp*) flint_malloc(sizeof(unity_zp) * (n_pow(2, k - 1) + 1));

    /* sets g_powers[0] = 1 */
    unity_zp_init(g_powers[0], f->p, f->exp, f->n);
    unity_zp_coeff_set_ui(g_powers[0], 0, 1);

    /* sets g_powers[1] = g */
    unity_zp_init(g_powers[1], f->p, f->exp, f->n);
    unity_zp_copy(g_powers[1], g);

    /* sets g_powers[i] = g^2 * g_powers[i - 1] */
    for (i = 2; i <= n_pow(2, k - 1); i++)
    {
        unity_zp_init(g_powers[i], f->p, f->exp, f->n);
        unity_zp_mul(g_powers[i], g_powers[i - 1], g_sqr);
    }

    unity_zp_set_zero(f);
    unity_zp_coeff_set_ui(f, 0, 1);
    i = fmpz_bits(pow) - 1;

    if (fmpz_equal_ui(pow, 1) == 1)
        unity_zp_copy(f, g);
    else {
    while (i >= 0)
    {
        if (digits[i] == 0)
        {
            unity_zp_sqr(temp, f);
            unity_zp_swap(temp, f);
            i--;
        }
        else
        {
            /* 
                finds length of chain; chain is length of
                longest bitstring less then k ending on 1.
            */
            j = i;
            chain = 0;
            while (j >= 0 && i - j < k) {
                if (digits[j] == 1) chain = i - j;
                j--;
            }
            /* f = f^(2^chain)) */
            for (h = 0; h <= chain; h++)
            {
                unity_zp_sqr(temp, f);
                unity_zp_swap(temp, f);
            }

            pow2k = 1;
            value = 0;

            /* 
                value = binary number (digits[i], ... , digits[i - chain])
                in decimal base
            */
            for (h = 0; h <= chain; h++)
            {
                value += digits[i - chain + h] * pow2k;
                pow2k *= 2;
            }

            /* f = f * g^value */
            unity_zp_mul(temp, f, g_powers[(value + 1) / 2]);
            unity_zp_swap(temp, f);

            /* increase i */
            i -= chain + 1;
        }
    }
    }

    for (i = 0; i <= n_pow(2, k - 1); i++)
        unity_zp_clear(g_powers[i]);
    flint_free(g_powers);
    flint_free(digits);

    fmpz_clear(power);
    unity_zp_clear(g_sqr);
    unity_zp_clear(temp);
}

