/*
    Copyright (C) 2009, 2010 William Hart
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"

void
_fmpz_vec_randtest(fmpz * f, flint_rand_t state, 
                   slong len, flint_bitcnt_t bits)
{
    slong i, sparseness;

    if (n_randint(state, 2))
    {
        for (i = 0; i < len; i++)
            fmpz_randtest(f + i, state, bits);
    }
    else
    {
        sparseness = 1 + n_randint(state, FLINT_MAX(2, len));

        for (i = 0; i < len; i++)
        {
            if (n_randint(state, sparseness))
                fmpz_zero(f + i);
            else
                fmpz_randtest(f + i, state, bits);
        }
    }
}

void
_fmpz_vec_randtest_unsigned(fmpz * f, flint_rand_t state, 
                            slong len, flint_bitcnt_t bits)
{
    slong i, sparseness;

    if (n_randint(state, 2))
    {
        for (i = 0; i < len; i++)
            fmpz_randtest_unsigned(f + i, state, bits);
    }
    else
    {
        sparseness = 1 + n_randint(state, FLINT_MAX(2, len));

        for (i = 0; i < len; i++)
        {
            if (n_randint(state, sparseness))
                fmpz_zero(f + i);
            else
                fmpz_randtest_unsigned(f + i, state, bits);
        }
    }
}
