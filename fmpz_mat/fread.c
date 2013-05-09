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

    Copyright (C) 2011  Andy Novocin

******************************************************************************/

#include <stdio.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mat.h"

int 
fmpz_mat_fread(FILE* file, fmpz_mat_t mat)
{
    len_t r, c, i, j;
    int byte_count;
    mpz_t t;

    /* first number in file should be row dimension */
    mpz_init(t);
    byte_count = mpz_inp_str(t, file, 10);
    if (byte_count == 0)
    {
        mpz_clear(t);
        return 0;
    }
    
    if (!mpz_fits_slong_p(t))
    {
        printf("Exception (fmpz_mat_fread). "
               "Number of rows does not fit into a len_t.\n");
        abort();
    }
    r = mpz_get_si(t);

    /* second number in file should be column dimension */
    byte_count = mpz_inp_str(t, file, 10);
    if (byte_count == 0)
    {
        mpz_clear(t);
        return 0;
    }
    
    if (!mpz_fits_slong_p(t))
    {
        printf("Exception (fmpz_mat_fread). "
               "Number of columns does not fit into a len_t.\n");
        abort();
    }
    c = mpz_get_si(t);
    mpz_clear(t);
    
    /* if the input is 0 by 0 then set the dimensions to r and c */
    if (mat->r == 0 && mat->c == 0)
    {
        fmpz_mat_clear(mat);
        fmpz_mat_init(mat,r,c);
    }
    else if (mat->r == 0 && mat->c == 0)
    {
        printf("Exception (fmpz_mat_fread). \n"
               "Dimensions are non-zero and do not match input dimensions.\n");
        abort();
    }

    for (i = 0; i < r; i++)
    {
        for (j = 0; j < c; j++)
        {
            if (!fmpz_fread(file, fmpz_mat_entry(mat, i, j)))
                return 0;
        }
    }

    /* a return value of 0 means a problem with 
       the file stream a value of 1 means success*/
    return 1;
}

