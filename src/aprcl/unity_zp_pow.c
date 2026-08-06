/*
    Copyright (C) 2015 Vladimir Glazachev

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_mod.h"
#include "aprcl.h"

/* -------------------------------------------------------------------------- */
/* from unity_zp_pow.c */

void
unity_zp_pow_fmpz(unity_zp f, const unity_zp g, const fmpz_t pow)
{
    slong i;
    unity_zp temp;

    unity_zp_init(temp, f->p, f->exp, fmpz_mod_ctx_modulus(f->ctx));

    unity_zp_set_zero(f);
    unity_zp_coeff_set_ui(f, 0, 1);

    i = fmpz_bits(pow) - 1;

    while (i >= 0)
    {
        unity_zp_sqr(temp, f);
        unity_zp_swap(f, temp);

        if (fmpz_tstbit(pow, i) == 1)
        {
            unity_zp_mul(temp, f, g);
            unity_zp_swap(f, temp);
        }

        i--;
    }

    unity_zp_clear(temp);
}

void
unity_zp_pow_ui(unity_zp f, const unity_zp g, ulong pow)
{
    fmpz_t p;
    fmpz_init_set_ui(p, pow);
    unity_zp_pow_fmpz(f, g, p);
    fmpz_clear(p);
}

/* -------------------------------------------------------------------------- */
/* from unity_zp_pow_2k.c */

void
unity_zp_pow_2k_fmpz(unity_zp f, const unity_zp g, const fmpz_t pow)
{
    ulong j, k, pow2k;
    slong i, e;
    fmpz_t digit;
    unity_zp temp;
    unity_zp *g_powers;

    fmpz_init(digit);
    unity_zp_init(temp, f->p, f->exp, fmpz_mod_ctx_modulus(f->ctx));

    /* g_sqr = g * g */
    unity_zp_sqr(temp, g);

    /* selects optimal k value for n */
    k = _unity_zp_pow_select_k(pow);
    /* selects e such that 2^(ek) < n < 2^((e + 1) * k) */
    e = (fmpz_bits(pow) - 1) / k;

    /*
        g_powers store odd powers of g up to 2^k - 1;
        g_powers[(i + 1) / 2] = g^i
    */
    pow2k = 1 << (k - 1);
    g_powers = (unity_zp*) flint_malloc(sizeof(unity_zp) * (pow2k + 1));

    /* sets g_powers[0] = 1 */
    unity_zp_init(g_powers[0], f->p, f->exp, fmpz_mod_ctx_modulus(f->ctx));
    unity_zp_coeff_set_ui(g_powers[0], 0, 1);

    /* sets g_powers[1] = g */
    unity_zp_init(g_powers[1], f->p, f->exp, fmpz_mod_ctx_modulus(f->ctx));
    unity_zp_copy(g_powers[1], g);

    /* sets g_powers[i] = g^2 * g_powers[i - 1] */
    for (i = 2; i <= pow2k; i++)
    {
        unity_zp_init(g_powers[i], f->p, f->exp, fmpz_mod_ctx_modulus(f->ctx));
        unity_zp_mul(g_powers[i], g_powers[i - 1], temp);
    }

    /* for all digits[i] */
    for (i = e; i >= 0; i--)
    {
        /*
            digit contains i-th digit of pow in k-ary base;
            k <= 11 so digit < 2^11 and fit into ulong
        */
        fmpz_fdiv_q_2exp(digit, pow, i * k);
        fmpz_fdiv_r_2exp(digit, digit, k);

        /* if digit == 0 set f = f^(2^k) */
        if (*digit == 0)
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

            /* digit = 2^t * b and b is odd */
            t = aprcl_p_power_in_q(*digit, 2);
            b = *digit / (1 << t);

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

    for (i = 0; i <= pow2k; i++)
        unity_zp_clear(g_powers[i]);
    flint_free(g_powers);

    fmpz_clear(digit);
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

/* -------------------------------------------------------------------------- */
/* from unity_zp_pow_select_k.c */

/*
    returns smallest integer k satisfies:
        log(n) < (k * (k + 1) * 2^(2 * k)) / (2^(k + 1) - k - 2) + 1
*/
ulong
_unity_zp_pow_select_k(const fmpz_t n)
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

/* -------------------------------------------------------------------------- */
/* from unity_zp_pow_sliding.c */

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
