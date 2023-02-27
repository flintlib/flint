/*
    Copyright (C) 2011, 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arith.h"

#if FLINT_BITS == 64
#define MAX_N_1LIMBS 25
#define MAX_N_2LIMBS 42
#define MAX_N_3LIMBS 58
#else
#define MAX_N_1LIMBS 15
#define MAX_N_2LIMBS 25
#define MAX_N_3LIMBS 34
#endif

static void
arith_bell_number_recursive(fmpz_t res, ulong n)
{
    mp_limb_t t[3 * MAX_N_3LIMBS];
    slong i, k;

    t[0] = 1;
    for (i = 1; i < FLINT_MIN(n, MAX_N_1LIMBS); i++)
    {
        t[i] = t[0];
        for (k = i; k >= 1; k--)
            t[k - 1] += t[k];
    }

    if (n <= MAX_N_1LIMBS)
    {
        fmpz_set_ui(res, t[0]);
        return;
    }

    for (k = i - 1; k >= 0; k--)
    {
        t[2 * k + 0] = t[k];
        t[2 * k + 1] = 0;
    }

    for ( ; i < FLINT_MIN(n, MAX_N_2LIMBS); i++)
    {
        t[2 * i + 0] = t[0];
        t[2 * i + 1] = t[1];

        for (k = i; k >= 1; k--)
            add_ssaaaa(t[2 * (k - 1) + 1], t[2 * (k - 1)],
                       t[2 * (k - 1) + 1], t[2 * (k - 1)],
                       t[2 * k + 1], t[2 * k]);
    }

    if (n <= MAX_N_2LIMBS)
    {
        fmpz_set_uiui(res, t[1], t[0]);
        return;
    }

    for (k = i - 1; k >= 0; k--)
    {
        t[3 * k + 2] = 0;
        t[3 * k + 1] = t[2 * k + 1];
        t[3 * k + 0] = t[2 * k];
    }

    for ( ; i < n; i++)
    {
        t[3 * i + 0] = t[0];
        t[3 * i + 1] = t[1];
        t[3 * i + 2] = t[2];

        for (k = i; k >= 1; k--)
            add_sssaaaaaa(t[3 * (k - 1) + 2], t[3 * (k - 1) + 1], t[3 * (k - 1)],
                          t[3 * (k - 1) + 2], t[3 * (k - 1) + 1], t[3 * (k - 1)],
                          t[3 * k + 2],       t[3 * k + 1],       t[3 * k]);
    }

    fmpz_set_ui_array(res, t, 3);
}

void
arith_bell_number(fmpz_t res, ulong n)
{
    if (n <= MAX_N_1LIMBS)
        fmpz_set_ui(res, bell_number_tab[n]);
    else if (n <= MAX_N_3LIMBS)
        arith_bell_number_recursive(res, n);
    else if (n <= 3400)
        arith_bell_number_dobinski(res, n);
    else
        arith_bell_number_multi_mod(res, n);
}
