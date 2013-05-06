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

    Copyright (C) 2010, 2011 Sebastian Pancratz

******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpq_poly.h"

char * _fmpq_poly_get_str_pretty(const fmpz *poly, const fmpz_t den, long len, 
                                 const char *var)
{
    long i;
    size_t j;
    size_t size;     /* Upper bound on the length          */
    size_t densize;  /* Size of the denominator in base 10 */
    size_t varsize;  /* Length of the variable name        */
    mpz_t z;         /* op->den (if this is not 1)         */
    mpq_t q;
    char *str;
    
    if (len == 0)  /* Zero polynomial */
    {
        str = flint_malloc(2);
        if (!str)
        {
            printf("Exception (fmpq_poly_get_str_pretty). malloc failed.\n");
            abort();
        }
        str[0] = '0';
        str[1] = '\0';
        return str;
    }
    if (len == 1)  /* Constant polynomials */
    {
        mpq_init(q);
        fmpz_get_mpz(mpq_numref(q), poly);
        fmpz_get_mpz(mpq_denref(q), den);
        mpq_canonicalize(q);
        str = mpq_get_str(NULL, 10, q);
        mpq_clear(q);
        return str;
    }
    if (len == 2)  /* Linear polynomials */
    {
        mpq_t a0, a1;
        size_t size0, size1;

        mpq_init(a0);
        mpq_init(a1);

        fmpz_get_mpz(mpq_numref(a0), poly);
        fmpz_get_mpz(mpq_denref(a0), den);
        mpq_canonicalize(a0);
        fmpz_get_mpz(mpq_numref(a1), poly + 1);
        fmpz_get_mpz(mpq_denref(a1), den);
        mpq_canonicalize(a1);

        size0 = mpz_sizeinbase(mpq_numref(a0), 10) 
              + mpz_sizeinbase(mpq_denref(a0), 10) + 1;
        size1 = mpz_sizeinbase(mpq_numref(a1), 10) 
              + mpz_sizeinbase(mpq_denref(a1), 10) + 1;
        size  = size0 + 1 + strlen(var) + 1 + size1 + 1;
        str   = flint_malloc(size);
        if (!str)
        {
            printf("Exception (fmpq_poly_get_str_pretty). malloc failed.\n");
            abort();
        }

        if (mpq_cmp_si(a1, 1, 1) == 0)
        {
            if (mpq_sgn(a0) == 0)
                gmp_sprintf(str, "%s", var);
            else if (mpq_sgn(a0) > 0)
                gmp_sprintf(str, "%s+%Qd", var, a0);
            else  /* mpq_sgn(a0) < 0 */
                gmp_sprintf(str, "%s%Qd", var, a0);
        }
        else if (mpq_cmp_si(a1, -1, 1) == 0)
        {
            if (mpq_sgn(a0) == 0)
                gmp_sprintf(str, "-%s", var);
            else if (mpq_sgn(a0) > 0)
                gmp_sprintf(str, "-%s+%Qd", var, a0);
            else  /* mpq_sgn(a0) < 0 */
                gmp_sprintf(str, "-%s%Qd", var, a0);
        }
        else
        {
            if (mpq_sgn(a0) == 0)
                gmp_sprintf(str, "%Qd*%s", a1, var);
            else if (mpq_sgn(a0) > 0)
                gmp_sprintf(str, "%Qd*%s+%Qd", a1, var, a0);
            else  /* mpq_sgn(a0) < 0 */
                gmp_sprintf(str, "%Qd*%s%Qd", a1, var, a0);
        }

        mpq_clear(a0);
        mpq_clear(a1);

        return str;
    }
    
    varsize = strlen(var);
    
    /* Copy the denominator into an mpz_t */
    mpz_init(z);
    if (*den == 1L)
    {
        densize = 0;
    }
    else
    {
        fmpz_get_mpz(z, den);
        densize = mpz_sizeinbase(z, 10);
    }
    
    /* Estimate the length */
    size = 0;
    for (i = 0; i < len; i++)
    {
        fmpz_get_mpz(z, poly + i);
        size += mpz_sizeinbase(z, 10);                   /* Numerator         */
        if (mpz_sgn(z) != 0) size += 1 + densize;        /* Denominator and / */
        size += 3;                                       /* Operator and ws   */
        size += 1 + varsize + 1;                         /* *, x and ^        */
        size += (size_t) ceil(log10((double) (i + 1)));  /* Exponent          */
    }
    
    mpq_init(q);
    str = flint_malloc(size);
    if (!str)
    {
        printf("Exception: malloc failed in fmpq_poly_to_string_pretty\n");
        mpz_clear(z);
        mpq_clear(q);
        abort();
    }
    j = 0;
    
    /* Print the leading term */
    fmpz_get_mpz(mpq_numref(q), poly + (len - 1));
    fmpz_get_mpz(mpq_denref(q), den);
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
    j += sprintf(str + j, "%li", len - 1);
    
    i = len - 1;
    while (i)
    {
        i--;
        
        if (fmpz_is_zero(poly + i))
            continue;
        
        fmpz_get_mpz(mpq_numref(q), poly + i);
        fmpz_get_mpz(mpq_denref(q), den);
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

char * fmpq_poly_get_str_pretty(const fmpq_poly_t poly, const char * var)
{
    return _fmpq_poly_get_str_pretty(poly->coeffs, poly->den, poly->length, 
                                     var);
}

