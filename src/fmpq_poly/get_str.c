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
#include <string.h>
#include <math.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpq_poly.h"

char * fmpq_poly_get_str(const fmpq_poly_t poly)
{
    slong i;
    size_t j;
    size_t len;     /* Upper bound on the length          */
    size_t denlen;  /* Size of the denominator in base 10 */
    mpz_t z;
    mpq_t q;
    char * str;
    
    if (poly->length == 0)
    {
        str = (char *) flint_malloc(2 * sizeof(char));
        str[0] = '0';
        str[1] = '\0';
        return str;
    }
    
    mpz_init(z);
    if (*poly->den == WORD(1))
    {
        denlen = 0;
    }
    else
    {
        fmpz_get_mpz(z, poly->den);
        denlen = mpz_sizeinbase(z, 10);
    }
    len = (size_t) ceil(log10((double) (poly->length + 1))) + (size_t) 2;
    for (i = 0; i < poly->length; i++)
    {
        fmpz_get_mpz(z, poly->coeffs + i);
        len += mpz_sizeinbase(z, 10) + (size_t) 1;
        if (mpz_sgn(z))
            len += denlen + (size_t) 2;
    }
    
    mpq_init(q);
    str = (char *) flint_malloc(len * sizeof(char));
    
    j = flint_sprintf(str, "%li", poly->length);
    str[j++] = ' ';
    for (i = 0; i < poly->length; i++)
    {
        str[j++] = ' ';
        fmpz_get_mpz(mpq_numref(q), poly->coeffs + i);
        fmpz_get_mpz(mpq_denref(q), poly->den);
        mpq_canonicalize(q);
        mpq_get_str(str + j, 10, q);
        j += strlen(str + j);
    }
    
    mpq_clear(q);
    mpz_clear(z);
    
    return str;
}

