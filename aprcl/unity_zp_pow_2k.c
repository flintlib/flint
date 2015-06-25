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

ulong
_unity_zp_pow_2k_find_k(const fmpz_t n)
{
    ulong bits;
    bits = fmpz_bits(n);

    if (bits <= 8)  return 1;
    if (bits <= 24) return 2;
    if (bits <= 69) return 3;
    if (bits <= 196) return 4;
    if (bits <= 538) return 5;
    if (bits <= 1433) return 6;
    if (bits <= 3714) return 7;
    if (bits <= 9399) return 8;
    if (bits <= 23290) return 9;
    if (bits <= 56651) return 10;
    return 11;
}

void
unity_zp_pow_2k_fmpz(unity_zp f, const unity_zp g, const fmpz_t pow)
{
    ulong j, k, pow2k;
    ulong *digits;
    slong i, e;
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
    /* selects e such that 2^(ek) < n < 2^((e + 1) * k) */
    e = (fmpz_bits(pow) - 1) / k;
    /* computes 2^k */
    pow2k = n_pow(2, k);

    /* 
        digits[i] = i-th digit of n in 2^k-ary base;
        most signifant digit have number 0
    */
    digits = (ulong*) flint_malloc(sizeof(ulong) * (e + 1));
    for (i = 0; i <= e; i++)
    {
        digits[i] = fmpz_tdiv_ui(power, pow2k);
        fmpz_tdiv_q_ui(power, power, pow2k);
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

    /* for all digits[i] */
    for (i = e; i >= 0; i--)
    {
        /* if digits[i] == 0 set f = f^(2^k) */
        if (digits[i] == 0)
        {
            for (j = 0; j < k; j++)
            {
                /* sets f = f^2 */
                unity_zp_sqr(temp, f);
                unity_zp_swap(temp, f);
            }
        }
        else
        {
            ulong t, b;

            /* digits[i] = 2^t * b and b is odd */
            t = p_power_in_q(digits[i], 2);
            b = digits[i] / n_pow(2, t);

            if (i == e) 
            {
                unity_zp_copy(f, g_powers[(b + 1) / 2]);
            }
            else
            {
                /* sets f = f^(2^(k - t)) */
                for (j = 0; j < k - t; j++)
                {
                    unity_zp_sqr(temp, f);
                    unity_zp_swap(temp, f);
                }

                /* sets f = f * g^b */
                unity_zp_mul(temp, f, g_powers[(b + 1) / 2]);
                unity_zp_swap(temp, f);
            }

            /* sets f = f^(2^t) */
            for (j = 0; j < t; j++)
            {
                unity_zp_sqr(temp, f);
                unity_zp_swap(temp, f);
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

void
unity_zp_pow_2k_ui(unity_zp f, const unity_zp g, ulong pow)
{
    fmpz_t p;
    fmpz_init_set_ui(p, pow);
    unity_zp_pow_2k_fmpz(f, g, p);
    fmpz_clear(p);
}

