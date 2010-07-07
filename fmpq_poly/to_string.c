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
    Copyright (C) 2010 William Hart
   
******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpq_poly.h"

ulong _fmpq_poly_decimal_digits(ulong n)
{
    ulong count;
    if (n == 0UL)
        return 1UL;
    count = 0UL;
    while (n > 0UL)
    {
        n /= 10UL;
        count++;
    }
    return count;
}

char * fmpq_poly_to_string(const fmpq_poly_t poly)
{
    ulong i, j;
    ulong len;     /* Upper bound on the length          */
    ulong denlen;  /* Size of the denominator in base 10 */
    mpz_t z;
    mpq_t q;
    
    char * str;
    
    if (poly->length == 0)
    {
        str = calloc(2, sizeof(char));
        str[0] = '0';
        str[1] = '\0';
        return str;
    }
    
    mpz_init(z);
    if (*poly->den == 1L)
    {
        denlen = 0;
    }
    else
    {
        fmpz_get_mpz(z, poly->den);
        denlen = mpz_sizeinbase(z, 10);
    }
    len = _fmpq_poly_decimal_digits(poly->length) + 2;
    for (i = 0; i < poly->length; i++)
    {
        fmpz_get_mpz(z, poly->coeffs + i);
        len += mpz_sizeinbase(z, 10) + 1;
        if (mpz_sgn(z))
            len += 2 + denlen;
    }
    
    mpq_init(q);
    str = calloc(len, sizeof(char));
    sprintf(str, "%lu", poly->length);
    for (j = 0; str[j] != '\0'; j++);
    str[j++] = ' ';
    for (i = 0; i < poly->length; i++)
    {
        str[j++] = ' ';
        fmpz_get_mpz(mpq_numref(q), poly->coeffs + i);
        fmpz_get_mpz(mpq_denref(q), poly->den);
        mpq_canonicalize(q);
        mpq_get_str(str + j, 10, q);
        for ( ; str[j] != '\0'; j++);
    }
    
    mpq_clear(q);
    mpz_clear(z);
    
    return str;
}

char * fmpq_poly_to_string_pretty(const fmpq_poly_t poly, const char * var)
{
    ulong i, j;
    ulong len;     /* Upper bound on the length          */
    ulong denlen;  /* Size of the denominator in base 10 */
    ulong varlen;  /* Length of the variable name        */
    mpz_t z;       /* op->den (if this is not 1)         */
    mpq_t q;
    char * str;
    
    if (poly->length == 0)  /* Zero polynomial */
    {
        str = calloc(2, sizeof(char));
        str[0] = '0';
        str[1] = '\0';
        return str;
    }
    if (poly->length == 1)  /* Constant polynomials */
    {
        mpq_init(q);
        fmpz_get_mpz(mpq_numref(q), poly->coeffs + (poly->length - 1UL));
        fmpz_get_mpz(mpq_denref(q), poly->den);
        mpq_canonicalize(q);
        str = mpq_get_str(NULL, 10, q);
        mpq_clear(q);
        return str;
    }
    
    varlen = strlen(var);
    
    /* Copy the denominator into an mpz_t */
    mpz_init(z);
    if (*poly->den == 1L)
    {
        denlen = 0;
    }
    else
    {
        fmpz_get_mpz(z, poly->den);
        denlen = mpz_sizeinbase(z, 10);
    }
    
    /* Estimate the length */
    len = 0;
    for (i = 0; i < poly->length; i++)
    {
        fmpz_get_mpz(z, poly->coeffs + i);
        len += mpz_sizeinbase(z, 10);            /* Numerator               */
        if (mpz_sgn(z) != 0) len += 1 + denlen;  /* Denominator and /       */
        len += 3;                                /* Operator and whitespace */
        len += 1 + varlen + 1;                   /* *, x and ^              */
        len += _fmpq_poly_decimal_digits(i);     /* Exponent                */
    }
    
    mpq_init(q);
    str = calloc(len, sizeof(char));
    j = 0;
    
    /* Print the leading term */
    fmpz_get_mpz(mpq_numref(q), poly->coeffs + (poly->length - 1UL));
    fmpz_get_mpz(mpq_denref(q), poly->den);
    mpq_canonicalize(q);
    
    if (mpq_cmp_si(q, 1, 1) != 0)
    {
        if (mpq_cmp_si(q, -1, 1) == 0)
            str[j++] = '-';
        else
        {
            mpq_get_str(str, 10, q);
            for ( ; str[j] != '\0'; j++);
            str[j++] = '*';
        }
    }
    sprintf(str + j, "%s", var);
    j += varlen;
    str[j++] = '^';
    sprintf(str + j, "%lu", poly->length - 1UL);
    for ( ; str[j] != '\0'; j++);
    
    i = poly->length - 1UL;
    while (i)
    {
        i--;
        
        if (fmpz_is_zero(poly->coeffs + i))
            continue;
        
        fmpz_get_mpz(mpq_numref(q), poly->coeffs + i);
        fmpz_get_mpz(mpq_denref(q), poly->den);
        mpq_canonicalize(q);
        
        str[j++] = ' ';
        if (mpq_sgn(q) < 0)
        {
            mpq_abs(q, q);
            str[j++] = '-';
        }
        else
            str[j++] = '+';
        str[j++] = ' ';
        
        mpq_get_str(str + j, 10, q);
        for ( ; str[j] != '\0'; j++);
        
        if (i > 0)
        {
            str[j++] = '*';
            sprintf(str + j, "%s", var);
            j += varlen;
            if (i > 1)
            {
                str[j++] = '^';
                sprintf(str + j, "%lu", i);
                for ( ; str[j] != '\0'; j++);
            }
        }
    }
    
    mpq_clear(q);
    mpz_clear(z);
    return str;
}

