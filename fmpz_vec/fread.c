/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"

int _fmpz_vec_fread(FILE * file, fmpz ** vec, len_t * len)
{
    int alloc, r;
    len_t i;
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
        printf("Exception (_fmpz_vec_fread). Length does not fit into a len_t.\n");
        abort();
    }
    if (alloc)
    {
        *len = mpz_get_si(t);
        *vec = _fmpz_vec_init(*len);
    }
    else
    {
        if (*len != mpz_get_si(t))
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

