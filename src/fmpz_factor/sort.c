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

// Tiny implem of quick sort, qsort is not usable because we also want to permute exp
// and qsort_r / qsort_s have portability issues
static void _fmpz_factor_sort_part(fmpz* p, ulong* exp, slong n)
{
    slong i, j;
    if (n <= 1) return;

    fmpz_swap(p + n / 2, p + n - 1);
    FLINT_SWAP(ulong, exp[n / 2], exp[n - 1]);
    for (i = 0, j = 0; i < n - 1; i++)
    {
        if (fmpz_cmp(p + i, p + n - 1) < 0)
        {
            fmpz_swap(p + i, p + j);
            FLINT_SWAP(ulong, exp[i], exp[j]);
            j++;
        }
    }
    fmpz_swap(p + j, p + n - 1);
    FLINT_SWAP(ulong, exp[j], exp[n - 1]);

    _fmpz_factor_sort_part(p, exp, j);
    _fmpz_factor_sort_part(p + j + 1, exp + j + 1, n - j - 1);
}

void
fmpz_factor_sort(fmpz_factor_t factor)
{
    _fmpz_factor_sort_part(factor->p, factor->exp, factor->num);
}
