/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "radix.h"

/* Functions to compute digit valuation. */

static ulong
_radix_val_exp1(ulong x, const void * FLINT_UNUSED(_radix))
{
    return 0;
}

static ulong
_radix_val_pow2(ulong x, const void * _radix)
{
    const radix_struct * radix = _radix;
    /* Todo: optimize this division */
    return flint_ctz(x) / radix->bval;
}

#define ODD_DIVISIBILITY_TEST(exp) \
    do { \
        ulong inv1 = radix->bpow_oddinv[exp].a; \
        ulong inv2 = radix->bpow_oddinv[exp].b; \
        if (n_divisible_odd_gm(x, inv1, inv2)) \
        { \
            x *= inv1; \
            val += exp; \
        } \
    } while (0)

static ulong
_radix_val_odd2(ulong x, const void * _radix)
{
    const radix_struct * radix = _radix;
    ulong val = 0;
    ODD_DIVISIBILITY_TEST(1);
    return val;
}

static ulong
_radix_val_odd4(ulong x, const void * _radix)
{
    const radix_struct * radix = _radix;
    ulong val = 0;
    ODD_DIVISIBILITY_TEST(2);
    ODD_DIVISIBILITY_TEST(1);
    return val;
}

static ulong
_radix_val_odd8(ulong x, const void * _radix)
{
    const radix_struct * radix = _radix;
    ulong val = 0;
    ODD_DIVISIBILITY_TEST(4);
    ODD_DIVISIBILITY_TEST(2);
    ODD_DIVISIBILITY_TEST(1);
    return val;
}

static ulong
_radix_val_odd16(ulong x, const void * _radix)
{
    const radix_struct * radix = _radix;
    ulong val = 0;
    ODD_DIVISIBILITY_TEST(8);
    ODD_DIVISIBILITY_TEST(4);
    ODD_DIVISIBILITY_TEST(2);
    ODD_DIVISIBILITY_TEST(1);
    return val;
}

static ulong
_radix_val_odd32(ulong x, const void * _radix)
{
    const radix_struct * radix = _radix;
    ulong val = 0;
    ODD_DIVISIBILITY_TEST(16);
    ODD_DIVISIBILITY_TEST(8);
    ODD_DIVISIBILITY_TEST(4);
    ODD_DIVISIBILITY_TEST(2);
    ODD_DIVISIBILITY_TEST(1);
    return val;
}

static ulong
_radix_val_odd64(ulong x, const void * _radix)
{
    const radix_struct * radix = _radix;
    ulong val = 0;
    ODD_DIVISIBILITY_TEST(32);
    ODD_DIVISIBILITY_TEST(16);
    ODD_DIVISIBILITY_TEST(8);
    ODD_DIVISIBILITY_TEST(4);
    ODD_DIVISIBILITY_TEST(2);
    ODD_DIVISIBILITY_TEST(1);
    return val;
}

static ulong
_radix_val_2odd2(ulong x, const void * _radix)
{
    const radix_struct * radix = _radix;
    unsigned int nz = flint_ctz(x);
    return FLINT_MIN(nz, _radix_val_odd2(x, radix));
}

static ulong
_radix_val_2odd4(ulong x, const void * _radix)
{
    const radix_struct * radix = _radix;
    unsigned int nz = flint_ctz(x);
    return FLINT_MIN(nz, _radix_val_odd4(x, radix));
}

static ulong
_radix_val_2odd8(ulong x, const void * _radix)
{
    const radix_struct * radix = _radix;
    unsigned int nz = flint_ctz(x);
    return FLINT_MIN(nz, _radix_val_odd8(x, radix));
}

static ulong
_radix_val_2odd16(ulong x, const void * _radix)
{
    const radix_struct * radix = _radix;
    unsigned int nz = flint_ctz(x);
    if (nz <= 1)
        return nz && n_divisible_odd_gm(x, radix->bpow_oddinv[1].a, radix->bpow_oddinv[1].b);
    return FLINT_MIN(nz, _radix_val_odd16(x, radix));
}

static ulong
_radix_val_2odd32(ulong x, const void * _radix)
{
    const radix_struct * radix = _radix;
    unsigned int nz = flint_ctz(x);
    if (nz <= 1)
        return nz && n_divisible_odd_gm(x, radix->bpow_oddinv[1].a, radix->bpow_oddinv[1].b);
    return FLINT_MIN(nz, _radix_val_odd32(x, radix));
}

/* This can be optimized, but there are probably few applications for
   radix b = 4*odd. */
static ulong
_radix_val_generic(ulong x, const void * _radix)
{
    const radix_struct * radix = _radix;
    unsigned int nz = flint_ctz(x);

    ulong val = 0;
    ulong inv1 = radix->bpow_oddinv[1].a;
    ulong inv2 = radix->bpow_oddinv[1].b;

    while (nz >= radix->bval && n_divisible_odd_gm(x, inv1, inv2))
    {
        x *= inv1;
        nz -= radix->bval;
        val++;
    }

    return val;
}

void radix_init(radix_t radix, ulong b, unsigned int exp)
{
    ulong B;
    ulong hi, lo;
    int i;

    if (b < 2 || exp >= FLINT_BITS)
        flint_throw(FLINT_ERROR, "radix_init: require b >= 2 and exp < FLINT_BITS");

    B = b;

    if (exp == 0)
    {
        for (i = 1; ; i++)
        {
            umul_ppmm(hi, lo, B, b);
            if (hi != 0)
            {
                exp = i;
                break;
            }
            else
            {
                B = lo;
            }
        }
    }
    else
    {
        for (i = 2; i <= exp; i++)
        {
            umul_ppmm(hi, B, B, b);
            if (hi != 0)
                flint_throw(FLINT_ERROR, "radix_init: require b^e < 2^FLINT_BITS");
        }
    }

    nmod_init(&radix->b, b);
    radix->exp = exp;
    nmod_init(&radix->B, B);

    /* todo: as a single malloc? */
    radix->bpow = flint_malloc(sizeof(ulong) * (exp + 1));
    radix->bpow_div = flint_malloc(sizeof(n_div_precomp_struct) * (exp + 1));

    /* todo: combine with earlier power computations */
    radix->bpow[0] = 1;
    radix->bpow[1] = b;
    for (i = 2; i <= exp; i++)
        radix->bpow[i] = radix->bpow[i - 1] * b;

    for (i = 0; i <= exp; i++)
        n_div_precomp_init(radix->bpow_div + i, radix->bpow[i]);

    {
        int prevbc = 0;
        for (i = 0; i <= exp; i++)
        {
            int jbc, bc = FLINT_BIT_COUNT(radix->bpow[i]);

            for (jbc = prevbc + 1; jbc <= bc; jbc++)
                radix->bits_to_digit_size[jbc - 1] = i;

            prevbc = bc;
        }
    }

    radix->bpow_oddinv = flint_malloc(sizeof(n_pair_struct) * (exp + 1));
    radix->bval = flint_ctz(b);
    radix->bpow_oddinv[0].a = 1;
    radix->bpow_oddinv[0].b = UWORD_MAX;

    ulong bodd_inv = n_binvert(b >> radix->bval);

    for (i = 1; i <= exp; i++)
    {
        ulong bodd = (radix->bpow[i] >> (radix->bval * i));

        radix->bpow_oddinv[i].a = radix->bpow_oddinv[i - 1].a * bodd_inv;
        radix->bpow_oddinv[i].b = UWORD_MAX / bodd;
    }

    if (exp == 1)
    {
        radix->val_func = _radix_val_exp1;
    }
    else if ((b & (b - 1)) == 0)
    {
        radix->val_func = _radix_val_pow2;
    }
    else if (radix->bval == 0)
    {
        if (exp == 2)       radix->val_func = _radix_val_odd2;
        else if (exp <= 4)  radix->val_func = _radix_val_odd4;
        else if (exp <= 8)  radix->val_func = _radix_val_odd8;
        else if (exp <= 16) radix->val_func = _radix_val_odd16;
        else if (exp <= 32) radix->val_func = _radix_val_odd32;
        else                radix->val_func = _radix_val_odd64;
    }
    else if (radix->bval == 1)  /* Mainly intended for base 10 */
    {
        if (exp == 2)       radix->val_func = _radix_val_2odd2;
        else if (exp <= 4)  radix->val_func = _radix_val_2odd4;
        else if (exp <= 8)  radix->val_func = _radix_val_2odd8;
        else if (exp <= 16) radix->val_func = _radix_val_2odd16;
        else                radix->val_func = _radix_val_2odd32;
    }
    else
    {
        radix->val_func = _radix_val_generic;
    }
}

void radix_clear(radix_t radix)
{
    flint_free(radix->bpow);
    flint_free(radix->bpow_div);
    flint_free(radix->bpow_oddinv);
}

