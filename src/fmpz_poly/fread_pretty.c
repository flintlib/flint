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
#include <ctype.h>
#include <gmp.h>

#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

static __inline__
int is_varsymbol0(char c)
{
    return isalpha((unsigned char) c);
}

static __inline__
int is_varsymbol1(char c)
{
    return isalnum((unsigned char) c) || (c == '_');
}

#define next_event()                                                      \
do {                                                                      \
    c = fgetc(file);                                                      \
    if (c == EOF && !feof(file))                                          \
        goto s_ioe;                                                       \
    ++r;                                                                  \
    if (i == N)                                                           \
    {                                                                     \
        buf = flint_realloc(buf, N = 2*N);                                \
    }                                                                     \
    buf[i++] = c;                                                         \
} while (0)

#define add_coeff()                                            \
do {                                                           \
    fmpz_set_mpz(f_coeff, z_coeff);                            \
    fmpz_poly_fit_length(poly, exp + 1);                       \
    fmpz_add(poly->coeffs + exp, poly->coeffs + exp, f_coeff); \
    if (poly->length < exp + 1)                                \
        poly->length = exp + 1;                                \
} while (0)

int fmpz_poly_fread_pretty(FILE *file, fmpz_poly_t poly, char **x)
{
    /*
        var     - variable name
        buf     - buffer of size N, at write position i, with c == buf[i-1]
        z_coeff - mpz_t for the coefficient, f_coeff is the fmpz_t version
        z_exp   - mpz_t for the exponent, exp is the slong version
        r       - return value
     */

    char *var = NULL;

    char c, *buf;
    int i, N;

    fmpz_t f_coeff;
    mpz_t z_coeff, z_exp;
    slong exp;

    int r = 0;

    fmpz_poly_zero(poly);
    if (poly->alloc)
        flint_mpn_zero((mp_ptr) poly->coeffs, poly->alloc);

    i = 0;
    N = 80;
    buf = flint_malloc(N);
    
    fmpz_init(f_coeff);
    mpz_init(z_coeff);
    mpz_init(z_exp);

    /* s_0 : */

    next_event();

    if (c == '-')
        goto s_1;
    if (isdigit((unsigned char) c))
        goto s_2;
    if (is_varsymbol0(c))
    {
        flint_mpz_set_si(z_coeff, 1);
        goto s_3;
    }

    goto s_parse_error;

  s_1 :

    next_event();

    if (isdigit((unsigned char) c))
        goto s_2;
    if (is_varsymbol0(c))
    {
        if (i == 1)
            flint_mpz_set_si(z_coeff, 1);
        else  /* i == 2 */
        {
            flint_mpz_set_si(z_coeff, -1);
            buf[0] = c;
            i = 1;
        }
        goto s_3;
    }

    goto s_parse_error;

  s_2 :

    next_event();

    if (isdigit((unsigned char) c))
        goto s_2;
    if (c == '*')
    {
        buf[i-1] = '\0';
        mpz_set_str(z_coeff, buf, 10);
        i = 0;
        goto s_4;
    }
    {
        buf[i-1] = '\0';
        mpz_set_str(z_coeff, buf, 10);
        exp = 0;
        add_coeff();
        if (c == '+' || c == '-')
        {
            i = 0;
            if (c == '-')
                buf[i++] = '-';
            goto s_1;
        }
        else
            goto s_end;
    }

  s_3 :

    next_event();

    if (is_varsymbol1(c))
        goto s_3;

    {
        buf[i-1] = '\0';
        if (var)
        {
            if (strcmp(buf, var))  /* Parse error */
                goto s_parse_error;
        }
        else
        {
            var = flint_malloc(i);
            strcpy(var, buf);
        }

        if (c == '^')
        {
            i = 0;
            goto s_5;
        }
        else if (c == '+' || c == '-')
        {
            exp = 1;
            add_coeff();

            i = 0;
            if (c == '-')
                buf[i++] = '-';
            goto s_1;
        }
        else
        {
            exp = 1;
            add_coeff();
            goto s_end;
        }
    }

  s_4 :

    next_event();

    if (is_varsymbol0(c))
        goto s_3;

    goto s_parse_error;

  s_5 :

    next_event();

    if (isdigit((unsigned char) c))
        goto s_6;

    goto s_parse_error;

  s_6 :

    next_event();

    if (isdigit((unsigned char) c))
        goto s_6;

    {
        buf[i-1] = '\0';
        mpz_set_str(z_exp, buf, 10);
        if (!mpz_fits_slong_p(z_exp))
        {
            goto s_parse_error;
        }
        exp = flint_mpz_get_si(z_exp);
        add_coeff();

        if (c == '+' || c == '-')
        {
            i = 0;
            if (c == '-')
                buf[i++] = '-';
            goto s_1;
        }
        else
            goto s_end;
    }

  s_parse_error :

    r = -2;
    goto s_end;

  s_ioe :

    r = -1;
    goto s_end;

  s_end :

    _fmpz_poly_normalise(poly);

    if (var)
        *x = var;
    else
    {
        *x = flint_malloc(1);
        **x = '\0';
    }

    fmpz_clear(f_coeff);
    mpz_clear(z_coeff);
    mpz_clear(z_exp);
    flint_free(buf);

    return r;
}

#undef next_event
#undef add_coeff

