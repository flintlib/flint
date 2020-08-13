/*
    Copyright (C) 2011  Andy Novocin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"

int 
fmpz_mat_fread(FILE* file, fmpz_mat_t mat)
{
    slong r, c, i, j;
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
        flint_printf("Exception (fmpz_mat_fread). "
               "Number of rows does not fit into a slong.\n");
        flint_abort();
    }
    r = flint_mpz_get_si(t);

    /* second number in file should be column dimension */
    byte_count = mpz_inp_str(t, file, 10);
    if (byte_count == 0)
    {
        mpz_clear(t);
        return 0;
    }
    
    if (!mpz_fits_slong_p(t))
    {
        flint_printf("Exception (fmpz_mat_fread). "
               "Number of columns does not fit into a slong.\n");
        flint_abort();
    }
    c = flint_mpz_get_si(t);
    mpz_clear(t);
    
    /* if the input is 0 by 0 then set the dimensions to r and c */
    if (mat->r == 0 && mat->c == 0)
    {
        fmpz_mat_clear(mat);
        fmpz_mat_init(mat,r,c);
    }
    else if (mat->r != r || mat->c != c)
    {
        flint_printf("Exception (fmpz_mat_fread). \n"
               "Dimensions are non-zero and do not match input dimensions.\n");
        flint_abort();
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

