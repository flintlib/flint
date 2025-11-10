/*
    Copyright (C) 2025 Mael Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_factor.h"

void
_fmpz_factor_gcd(fmpz_factor_t factor, const fmpz_factor_t a, const fmpz_factor_t b)
{
    slong i, j;
    ulong num_len;
    int cmp;

    num_len = FLINT_MIN(a->num, b->num);
    _fmpz_factor_fit_length(factor, num_len);
    factor->sign = 1;
    factor->num = 0;

    i = 0;
    j = 0;
    factor->num = 0;

    while ((i < a->num) && (j < b->num))
    {
        cmp = fmpz_cmp(a->p + i, b->p + j);
        if (cmp < 0) // a[i] < b[i]
        {
            i++;
        }
        else if (cmp > 0) // b[i] < a[i]
        {
            j++;
        }
        else // b[i] = a[i]
        {
            fmpz_set(factor->p + factor->num, a->p + i);
            factor->exp[factor->num] = FLINT_MIN(a->exp[i], b->exp[j]);
            factor->num++; 
            i++;
            j++;
        }
    }
}

void fmpz_factor_gcd(fmpz_factor_t factor, fmpz_factor_t a, fmpz_factor_t b)
{
    fmpz_factor_sort(a);
    fmpz_factor_sort(b);

    _fmpz_factor_gcd(factor, a, b);
}