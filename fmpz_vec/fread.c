/*
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"

int _fmpz_vec_fread(FILE * file, fmpz ** vec, slong * len)
{
    int alloc, r;
    slong i;
    mpz_t t;

    alloc = (*vec == NULL);

    mpz_init(t);
    r = mpz_inp_str(t, file, 10);
    if (r == 0)
    {
        if (alloc)
            *len = 0;
        mpz_clear(t);
        return 0;
    }
    if (!mpz_fits_slong_p(t))
    {
        flint_printf("Exception (_fmpz_vec_fread). Length does not fit into a slong.\n");
        flint_abort();
    }
    if (alloc)
    {
        *len = flint_mpz_get_si(t);
        *vec = _fmpz_vec_init(*len);
    }
    else
    {
        if (*len != flint_mpz_get_si(t))
        {
            mpz_clear(t);
            return 0;
        }
    }
    mpz_clear(t);

    for (i = 0; i < *len; i++)
    {
        r = fmpz_fread(file, (*vec) + i);
        if (r <= 0)
        {
            if (alloc)
            {
                _fmpz_vec_clear(*vec, *len);
                *vec = NULL;
                *len = 0;
            }
            return r;
        }
    }

    return 1;
}

