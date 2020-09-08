/*
    Copyright (C) 2015 Vladimir Glazachev

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "aprcl.h"

void
unity_zp_pow_sliding_fmpz(unity_zp f, unity_zp g, const fmpz_t pow)
{
    ulong h, k, value;
    slong i, j;
    unity_zp temp;
    unity_zp *g_powers;

    fmpz_t * t;
    t = (fmpz_t*) flint_malloc(sizeof(fmpz_t) * SQUARING_SPACE);
    for (i = 0; i < SQUARING_SPACE; i++)
        fmpz_init(t[i]);

    unity_zp_init(temp, f->p, f->exp, fmpz_mod_ctx_modulus(f->ctx));

    /* reduce g by cyclopoly */
    _unity_zp_reduce_cyclotomic(g);

    /* temp = g * g */
    unity_zp_sqr_inplace(temp, g, t);

    /* selects optimal k value for n */
    k = _unity_zp_pow_select_k(pow);

    /* 
        g_powers store odd powers of g up to 2^k - 1;
        g_powers[(i + 1) / 2] = g^i
    */
    g_powers = (unity_zp*) flint_malloc(sizeof(unity_zp) * (n_pow(2, k - 1) + 1));

    /* sets g_powers[0] = 1 */
    unity_zp_init(g_powers[0], f->p, f->exp, fmpz_mod_ctx_modulus(f->ctx));
    unity_zp_coeff_set_ui(g_powers[0], 0, 1);

    /* sets g_powers[1] = g */
    unity_zp_init(g_powers[1], f->p, f->exp, fmpz_mod_ctx_modulus(f->ctx));
    unity_zp_copy(g_powers[1], g);

    /* sets g_powers[i] = g^2 * g_powers[i - 1] */
    for (i = 2; i <= n_pow(2, k - 1); i++)
    {
        unity_zp_init(g_powers[i], f->p, f->exp, fmpz_mod_ctx_modulus(f->ctx));
        unity_zp_mul_inplace(g_powers[i], g_powers[i - 1], temp, t);
    }

    unity_zp_set_zero(f);
    unity_zp_coeff_set_ui(f, 0, 1);
    i = fmpz_bits(pow) - 1;

    /* working with pow = (e_l, e_{l-1}, ... , e_0) in 2 base */
    while (i >= 0)
    {
        if (fmpz_tstbit(pow, i) == 0)
        {
            unity_zp_sqr_inplace(temp, f, t);
            unity_zp_swap(temp, f);
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
                unity_zp_sqr_inplace(temp, f, t);
                unity_zp_swap(temp, f);
            }

            /* 
                value = binary number (e_i, ... , e_j) in decimal base
            */
            value = 0;
            for (h = 0; h < i - j + 1; h++)
                value += fmpz_tstbit(pow, j + h) << h;

            /* f = f * g^value */
            unity_zp_mul_inplace(temp, f, g_powers[(value + 1) / 2], t);
            unity_zp_swap(temp, f);

            /* increase i */
            i = j - 1;
        }
    }

    for (i = 0; i < SQUARING_SPACE; i++)
        fmpz_clear(t[i]);
    flint_free(t);

    for (i = 0; i <= n_pow(2, k - 1); i++)
        unity_zp_clear(g_powers[i]);
    flint_free(g_powers);

    unity_zp_clear(temp);
}

