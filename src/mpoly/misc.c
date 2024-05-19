/*
    Copyright (C) 2016-2017 William Hart
    Copyright (C) 2017-2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"

/* choose m so that (m + 1)/(n - m) ~= la/lb, i.e. m = (n*la - lb)/(la + lb) */
slong mpoly_divide_threads(slong n, double la, double lb)
{
    double m_double = (n*la - lb)/(la + lb);
    slong m = m_double + (2*m_double > n ? -0.5 : 0.5);

    /* input must satisfy */
    FLINT_ASSERT(n > 0);

    if (m <= 0)
        m = 0;

    if (m >= n - 1)
        m = n - 1;

    /* output must satisfy */
    FLINT_ASSERT(m >= 0);
    FLINT_ASSERT(m < n);
    return m;
}

void mpoly_monomial_msub_ui_array(
        ulong * exp1, const ulong * exp2,
        const ulong * scalar, slong scalar_limbs,
        const ulong * exp3, slong N)
{
    slong i;

    /* TODO: Check if exp1 == exp2? */
    for (i = 0; i < N; i++)
        exp1[i] = exp2[i];

    FLINT_ASSERT(scalar_limbs <= N);

    for (i = 0; i < scalar_limbs; i++)
    {
        FLINT_ASSERT(N > i);
        __gmpn_submul_1(exp1 + i, exp3, N - i, scalar[i]);
    }
}

void mpoly_monomial_madd_ui_array(
        ulong * exp1, const ulong * exp2,
        const ulong * scalar, slong scalar_limbs,
        const ulong * exp3, slong N)
{
    slong i;

    /* TODO: Check if exp1 == exp2? */
    for (i = 0; i < N; i++)
        exp1[i] = exp2[i];

    FLINT_ASSERT(scalar_limbs <= N);

    for (i = 0; i < scalar_limbs; i++)
        __gmpn_addmul_1(exp1 + i, exp3, N - i, scalar[i]);
}

void mpoly_monomial_max(ulong * exp1, const ulong * exp2, const ulong * exp3,
        flint_bitcnt_t bits, slong N, ulong mask)
{
    ulong s, m;
    slong i;
    for (i = 0; i < N; i++)
    {
        s = mask + exp2[i] - exp3[i];
        m = mask & s;
        m = m - (m >> (bits - 1));
        exp1[i] = exp3[i] + (s & m);
    }
}

void mpoly_monomial_min(ulong * exp1, const ulong * exp2, const ulong * exp3,
        flint_bitcnt_t bits, slong N, ulong mask)
{
    ulong s, m;
    slong i;
    for (i = 0; i < N; i++)
    {
        s = mask + exp2[i] - exp3[i];
        m = mask & s;
        m = m - (m >> (bits - 1));
        exp1[i] = exp2[i] - (s & m);
    }
}

void mpoly_monomial_max_mp(ulong * exp1, const ulong * exp2, const ulong * exp3,
        flint_bitcnt_t bits, slong N)
{
    slong i;
    flint_bitcnt_t j;
    for (i = 0; i < N; i += bits/FLINT_BITS)
    {
        const ulong * t = exp2;
        for (j = bits/FLINT_BITS - 1; (slong) j >= 0; j--)
        {
            if (exp3[i + j] != exp2[i + j])
            {
                if (exp3[i + j] > exp2[i + j])
                    t = exp3;
                break;
            }
        }
        for (j = 0; j < bits/FLINT_BITS; j++)
        {
            exp1[i + j] = t[i + j];
        }
    }
}

void mpoly_monomial_min_mp(ulong * exp1, const ulong * exp2, const ulong * exp3,
        flint_bitcnt_t bits, slong N)
{
    slong i;
    flint_bitcnt_t j;
    for (i = 0; i < N; i += bits/FLINT_BITS)
    {
        const ulong * t = exp2;
        for (j = bits/FLINT_BITS - 1; (slong) j >= 0; j--)
        {
            if (exp3[i + j] != exp2[i + j])
            {
                if (exp3[i + j] < exp2[i + j])
                    t = exp3;
                break;
            }
        }
        for (j = 0; j < bits/FLINT_BITS; j++)
        {
            exp1[i + j] = t[i + j];
        }
    }
}

void mpoly_main_variable_terms1(slong * i1, slong * n1, const ulong * exp1,
        slong l1, slong len1, slong k, slong FLINT_UNUSED(num), slong bits)
{
    slong i, j = 0;
    slong shift = bits*(k - 1);

    i1[0] = 0;
    for (i = 0; i < l1 - 1; i++)
    {
        while (j < len1 && (l1 - i - 1) == (slong) (exp1[j] >> shift))
            j++;

        i1[i + 1] = j;
        n1[i] = j - i1[i];
    }
    n1[l1 - 1] = len1 - j;
}
