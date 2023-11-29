/*
    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <gmpcompat.h>
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

/* printing *******************************************************************/

int
_fmpz_poly_fprint_pretty(FILE * file, const fmpz * poly, slong len, const char * x)
{
    int r;
    slong i;

    FMPZ_VEC_NORM(poly, len);

    if (len == 0)
    {
        r = fputc('0', file);
        r = (r != EOF) ? 1 : EOF;
        return r;
    }
    else if (len == 1)
    {
        r = fmpz_fprint(file, poly + 0);
        return r;
    }
    else if (len == 2)
    {
        if (*(poly + 1) == WORD(1))
        {
            r = flint_fprintf(file, "%s", x);
        }
        else if (*(poly + 1) == WORD(-1))
        {
            r = flint_fprintf(file, "-%s", x);
        }
        else
        {
            r = fmpz_fprint(file, poly + 1);
            if (r > 0)
                r = flint_fprintf(file, "*%s", x);
        }

        if (r > 0)
        {
            if (fmpz_sgn(poly + 0) > 0)
            {
                r = flint_fprintf(file, "+");
                if (r > 0)
                    r = fmpz_fprint(file, poly + 0);
            }
            else if (fmpz_sgn(poly + 0) < 0)
            {
                r = fmpz_fprint(file, poly + 0);
            }
        }
        return r;
    }

    i = len - 1;  /* i >= 2 */
    r = 1;
    {
        if (*(poly + i) == 1)
           r = flint_fprintf(file, "%s^%wd", x, i);
        else if (*(poly + i) == -1)
           r = flint_fprintf(file, "-%s^%wd", x, i);
        else
        {
           r = fmpz_fprint(file, poly + i);
           if (r > 0)
              r = flint_fprintf(file, "*%s^%wd", x, i);
        }
        --i;
    }

    for (; (r > 0) && (i > 1); --i)
    {
        if (*(poly + i) == 0)
            continue;

        if (*(poly + i) == 1)
            r = flint_fprintf(file, "+%s^%wd", x, i);
        else if (*(poly + i) == -1)
            r = flint_fprintf(file, "-%s^%wd", x, i);
        else
        {
            if (fmpz_sgn(poly + i) > 0)
            {
                r = fputc('+', file);
                r = (r != EOF) ? 1 : EOF;
            }
            if (r > 0)
                r = fmpz_fprint(file, poly + i);
            if (r > 0)
                r = flint_fprintf(file, "*%s^%wd", x, i);
        }
    }

    if ((r > 0) && *(poly + 1))
    {
        if (*(poly + 1) == 1)
        {
            r = fputc('+', file);
            r = (r != EOF) ? 1 : EOF;
            if (r > 0)
            {
                r = fputs(x, file);
                r = (r >= 0) ? 1 : -1;
            }
        }
        else if (*(poly + 1) == -1)
        {
            r = fputc('-', file);
            r = (r != EOF) ? 1 : EOF;
            if (r > 0)
            {
                r = fputs(x, file);
                r = (r >= 0) ? 1 : -1;
            }
        }
        else
        {
            if (fmpz_sgn(poly + 1) > 0)
            {
                r = fputc('+', file);
                r = (r != EOF) ? 1 : EOF;
            }
            if (r > 0)
                r = fmpz_fprint(file, poly + 1);
            if (r > 0)
            {
                r = fputc('*', file);
                r = (r != EOF) ? 1 : EOF;
            }
            if (r > 0)
            {
                r = fputs(x, file);
                r = (r >= 0) ? 1 : -1;
            }
        }
    }
    if ((r > 0) && *(poly))
    {
        if (fmpz_sgn(poly) > 0)
        {
            r = fputc('+', file);
            r = (r != EOF) ? 1 : EOF;
        }
        if (r > 0)
            r = fmpz_fprint(file, poly);
    }

    return r;
}

int _fmpz_poly_fprint(FILE * file, const fmpz * poly, slong len) { return _fmpz_vec_fprint(file, poly, len); }
int fmpz_poly_fprint(FILE * file, const fmpz_poly_t poly) { return _fmpz_vec_fprint(file, poly->coeffs, poly->length); }

int fmpz_poly_fprint_pretty(FILE * file, const fmpz_poly_t poly, const char * x) { return _fmpz_poly_fprint_pretty(file, poly->coeffs, poly->length, x); }

int _fmpz_poly_print(const fmpz * poly, slong n) { return _fmpz_poly_fprint(stdout, poly, n); }
int fmpz_poly_print(const fmpz_poly_t poly) { return fmpz_poly_fprint(stdout, poly); }

int _fmpz_poly_print_pretty(const fmpz * poly, slong len, const char * x) { return _fmpz_poly_fprint_pretty(stdout, poly, len, x); }
int fmpz_poly_print_pretty(const fmpz_poly_t poly, const char * x) { return fmpz_poly_fprint_pretty(stdout, poly, x); }

/* reading ********************************************************************/

int fmpz_poly_fread(FILE * file, fmpz_poly_t poly)
{
    int r;
    slong i, len;
    mpz_t t;

    mpz_init(t);
    r = mpz_inp_str(t, file, 10);
    if (r == 0)
    {
        mpz_clear(t);
        return 0;
    }
    if (!mpz_fits_slong_p(t))
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_poly_fread). Length does not fit into a slong.\n");
    }
    len = flint_mpz_get_si(t);
    mpz_clear(t);

    fmpz_poly_fit_length(poly, len);

    for (i = 0; i < len; i++)
    {
        r = fmpz_fread(file, poly->coeffs + i);
        if (r <= 0)
            return r;
    }

    _fmpz_poly_set_length(poly, len);
    _fmpz_poly_normalise(poly);

    return 1;
}

static inline
int is_varsymbol0(char c)
{
    return isalpha((unsigned char) c);
}

static inline
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

int fmpz_poly_read(fmpz_poly_t poly) { return fmpz_poly_fread(stdin, poly); }
int fmpz_poly_read_pretty(fmpz_poly_t poly, char **x) { return fmpz_poly_fread_pretty(stdin, poly, x); }

/* debugging ******************************************************************/

void fmpz_poly_debug(const fmpz_poly_t poly)
{
    flint_printf("(alloc = %wd, length = %wd, vec = ", poly->alloc, poly->length);
    if (poly->coeffs)
    {
        flint_printf("{");
        _fmpz_vec_print(poly->coeffs, poly->alloc);
        flint_printf("}");
    }
    else
    {
        flint_printf("NULL");
    }
    flint_printf(")");
    fflush(stdout);
}
