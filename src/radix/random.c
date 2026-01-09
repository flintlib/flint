/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "radix.h"

void
radix_init_randtest(radix_t radix, flint_rand_t state)
{
    switch (n_randint(state, 10))
    {
        case 0: radix_init(radix, 10, 1); break;
        case 1: radix_init(radix, 10, 2); break;
        case 2: radix_init(radix, 10, 3); break;
        case 3: radix_init(radix, 2 + n_randtest(state) % (UWORD_MAX - 1), 0); break;
        case 4: radix_init(radix, 2 + n_randtest(state) % (UWORD_MAX - 1), 1); break;
        case 5: radix_init(radix, UWORD_MAX, 1); break;
        case 6: radix_init(radix, 2, 1); break;
        case 7: radix_init(radix, 3, 1); break;
        default: radix_init(radix, 10, 0);
    }
}

void
radix_rand_limbs(nn_ptr res, flint_rand_t state, slong n, const radix_t radix)
{
    slong i;
    for (i = 0; i < n; i++)
        res[i] = n_urandint(state, LIMB_RADIX(radix));
}

void
radix_rand_digits(nn_ptr res, flint_rand_t state, slong n, const radix_t radix)
{
    slong nr = (n + radix->exp - 1) / radix->exp;

    radix_rand_limbs(res, state, nr, radix);
    if (n % radix->exp != 0)
        res[nr - 1] %= n_pow(DIGIT_RADIX(radix), n % radix->exp); /* todo: precomp */
}

void
radix_randtest_limbs(nn_ptr res, flint_rand_t state, slong n, const radix_t radix)
{
    slong start, run, i;
    ulong B = LIMB_RADIX(radix);

    start = 0;

    while (start < n)
    {
        run = 1 + n_randint(state, n - start);

        switch (n_randint(state, 3))
        {
            case 0:
                for (i = 0; i < run; i++)
                    res[start + i] = 0;
                break;
            case 1:
                for (i = 0; i < run; i++)
                    res[start + i] = B - 1;
                break;
            default:
                for (i = 0; i < run; i++)
                    res[start + i] = n_randint(state, B);
        }

        start += run;
    }
}

void
radix_randtest_digits(nn_ptr res, flint_rand_t state, slong n, const radix_t radix)
{
    slong nr = (n + radix->exp - 1) / radix->exp;

    radix_randtest_limbs(res, state, nr, radix);
    if (n % radix->exp != 0)
        res[nr - 1] %= n_pow(DIGIT_RADIX(radix), n % radix->exp); /* todo: precomp */
}

