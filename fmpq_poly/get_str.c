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
#include <string.h>
#include <math.h>
#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpq_poly.h"

char * fmpq_poly_get_str(const fmpq_poly_t poly)
{
    long i;
    size_t j;
    size_t len;     /* Upper bound on the length          */
    size_t denlen;  /* Size of the denominator in base 10 */
    mpz_t z;
    mpq_t q;
    char * str;
    
    if (poly->length == 0)
    {
        str = (char *) malloc(2 * sizeof(char));
        if (str == NULL)
        {
            printf("Exception: malloc failed in fmpq_poly_to_string\n");
            abort();
        }
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
    len = (size_t) ceil(log10((double) (poly->length + 1))) + (size_t) 2;
    for (i = 0; i < poly->length; i++)
    {
        fmpz_get_mpz(z, poly->coeffs + i);
        len += mpz_sizeinbase(z, 10) + (size_t) 1;
        if (mpz_sgn(z))
            len += denlen + (size_t) 2;
    }
    
    mpq_init(q);
    str = (char *) malloc(len * sizeof(char));
    if (str == NULL)
    {
        printf("Exception: malloc failed in fmpq_poly_to_string\n");
        mpz_clear(z);
        mpq_clear(q);
        abort();
    }
    j = sprintf(str, "%li", poly->length);
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

char * fmpq_poly_get_str_pretty(const fmpq_poly_t poly, const char * var)
{
    long i;
    size_t j;
    size_t len;     /* Upper bound on the length          */
    size_t denlen;  /* Size of the denominator in base 10 */
    size_t varlen;  /* Length of the variable name        */
    mpz_t z;        /* op->den (if this is not 1)         */
    mpq_t q;
    char * str;
    
    if (poly->length == 0)  /* Zero polynomial */
    {
        str = (char *) malloc(2 * sizeof(char));
        if (str == NULL)
        {
            printf("Exception: malloc failed in fmpq_poly_to_string_pretty\n");
            abort();
        }
        str[0] = '0';
        str[1] = '\0';
        return str;
    }
    if (poly->length == 1)  /* Constant polynomials */
    {
        mpq_init(q);
        fmpz_get_mpz(mpq_numref(q), poly->coeffs + (poly->length - 1));
        fmpz_get_mpz(mpq_denref(q), poly->den);
        mpq_canonicalize(q);
        str = mpq_get_str(NULL, 10, q);
        mpq_clear(q);
        return str;
    }
    if (poly->length == 2)  /* Degree 1 polynomials */
    {
        mpq_t a0, a1;
        size_t len0, len1;

        mpq_init(a0);
        mpq_init(a1);

        fmpz_get_mpz(mpq_numref(a0), poly->coeffs);
        fmpz_get_mpz(mpq_denref(a0), poly->den);
        mpq_canonicalize(a0);
        fmpz_get_mpz(mpq_numref(a1), poly->coeffs + 1);
        fmpz_get_mpz(mpq_denref(a1), poly->den);
        mpq_canonicalize(a1);

        len0 = mpz_sizeinbase(mpq_numref(a0), 10) 
             + mpz_sizeinbase(mpq_denref(a0), 10) + 1;
        len1 = mpz_sizeinbase(mpq_numref(a1), 10) 
             + mpz_sizeinbase(mpq_denref(a1), 10) + 1;
        len  = len0 + 1 + strlen(var) + 1 + len1 + 1;
        str  = malloc(len);
        if (!str)
        {
            printf("Exception: malloc failed in fmpq_poly_to_string_pretty\n");
            abort();
        }

        if (mpq_sgn(a0) == 0)
        {
            gmp_sprintf(str, "%Qd*%s", a1, var);
        }
        else if (mpq_sgn(a0) > 0)
        {
            gmp_sprintf(str, "%Qd*%s+%Qd", a1, var, a0);
        }
        else  /* mpq_sgn(a0) < 0 */
        {
            gmp_sprintf(str, "%Qd*%s%Qd", a1, var, a0);
        }

        mpq_clear(a0);
        mpq_clear(a1);

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
        len += mpz_sizeinbase(z, 10);                   /* Numerator         */
        if (mpz_sgn(z) != 0) len += 1 + denlen;         /* Denominator and / */
        len += 3;                                       /* Operator and ws   */
        len += 1 + varlen + 1;                          /* *, x and ^        */
        len += (size_t) ceil(log10((double) (i + 1)));  /* Exponent          */
    }
    
    mpq_init(q);
    str = (char *) malloc(len * sizeof(char));
    if (str == NULL)
    {
        printf("Exception: malloc failed in fmpq_poly_to_string_pretty\n");
        mpz_clear(z);
        mpq_clear(q);
        abort();
    }
    j = 0;
    
    /* Print the leading term */
    fmpz_get_mpz(mpq_numref(q), poly->coeffs + (poly->length - 1));
    fmpz_get_mpz(mpq_denref(q), poly->den);
    mpq_canonicalize(q);
    
    if (mpq_cmp_si(q, 1, 1) != 0)
    {
        if (mpq_cmp_si(q, -1, 1) == 0)
            str[j++] = '-';
        else
        {
            mpq_get_str(str, 10, q);
            j += strlen(str + j);
            str[j++] = '*';
        }
    }
    j += sprintf(str + j, "%s", var);
    str[j++] = '^';
    j += sprintf(str + j, "%li", poly->length - 1);
    
    i = poly->length - 1;
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
        j += strlen(str + j);
        
        if (i > 0)
        {
            str[j++] = '*';
            j += sprintf(str + j, "%s", var);
            if (i > 1)
            {
                str[j++] = '^';
                j += sprintf(str + j, "%li", i);
            }
        }
    }
    
    mpq_clear(q);
    mpz_clear(z);
    return str;
}

