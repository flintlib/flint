/*
    Copyright (C) 2010, 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "gmpcompat.h"
#include "fmpq.h"
#include "fmpq_poly.h"

/* printing *******************************************************************/

/*
    Recall the return value conventions for fputc (of type int)

    ``If there are no errors, the same character that has been written is
    returned.  If an error occurs, EOF is returned and the error indicator
    is set''

    where the EOF macro expands to a negative int, and flint_fprintf (of type int)

    ``On success, the total number of characters written is returned.
    On failure, a negative number is returned.''
 */

/* TODO: Remove me! */
static void
___fmpq_poly_set_array_mpq(fmpz * poly, fmpz_t den, const mpq_t * a, slong n)
{
    slong i;
    mpz_t d, t;

    flint_mpz_init_set_ui(d, 1);
    mpz_init(t);
    for (i = 0; i < n; i++)
    {
        mpz_lcm(d, d, mpq_denref(a[i]));
    }

    for (i = 0; i < n; i++)
    {
        __mpz_struct *ptr = _fmpz_promote(poly + i);

        mpz_divexact(t, d, mpq_denref(a[i]));
        mpz_mul(ptr, mpq_numref(a[i]), t);
        _fmpz_demote_val(poly + i);
    }

    fmpz_set_mpz(den, d);
    mpz_clear(d);
    mpz_clear(t);
}

/* TODO: Remove me! */
static void __fmpq_poly_set_array_mpq(fmpq_poly_t poly, const mpq_t * a, slong n)
{
    if (n == 0)
    {
        fmpq_poly_zero(poly);
    }
    else
    {
        fmpq_poly_fit_length(poly, n);
        ___fmpq_poly_set_array_mpq(poly->coeffs, poly->den, a, n);
        _fmpq_poly_set_length(poly, n);
        _fmpq_poly_normalise(poly);
    }
}

int
_fmpq_poly_fprint(FILE * file, const fmpz * poly, const fmpz_t den, slong len)
{
    int r;
    slong i;
    fmpz_t n, d, g;

    fmpz_init(n);
    fmpz_init(d);
    fmpz_init(g);

    r = flint_fprintf(file, "%wd", len);
    if ((len > 0) && (r > 0))
    {
        r = fputc(' ', file);
        for (i = 0; (i < len) && (r > 0); i++)
        {
            r = fputc(' ', file);
            if (r > 0)
            {
                fmpz_gcd(g, poly + i, den);
                fmpz_divexact(n, poly + i, g);
                fmpz_divexact(d, den, g);
                if (*d == WORD(1))
                    r = fmpz_fprint(file, n);
                else
                {
                    r = fmpz_fprint(file, n);
                    if (r > 0)
                        r = fputc('/', file);
                    if (r > 0)
                        r = fmpz_fprint(file, d);
                }
            }
        }
    }

    fmpz_clear(n);
    fmpz_clear(d);
    fmpz_clear(g);

    return r;
}

/*
    Macro wrapping _fmpq_fprint(file, x, y), ensuring that the printed
    rational is in lowest terms.  Assumes that y > 0.
 */

#define __fmpq_fprint(x,y)         \
do {                               \
    fmpz_gcd(g, x, y);             \
    if (fmpz_is_one(g))            \
    {                              \
        _fmpq_fprint(file, x, y);  \
    }                              \
    else                           \
    {                              \
        fmpz_divexact(n, x, g);    \
        fmpz_divexact(d, y, g);    \
        _fmpq_fprint(file, n, d);  \
    }                              \
} while (0)

/* checks if x/y == 1, where (x, y) need not be in lowest terms */
#define __fmpq_is_one(x,y) fmpz_equal((x), (y))

/* checks if x/y == +/- 1, where (x, y) need not be in lowest terms */
#define __fmpq_is_pm1(x,y) (fmpz_cmpabs((x),(y)) == 0)

int _fmpq_poly_fprint_pretty(FILE * file,
                             const fmpz *poly, const fmpz_t den, slong len,
                             const char * x)
{
    fmpz_t n, d, g;

    fmpz_init(n);
    fmpz_init(d);
    fmpz_init(g);

    if (len == 0)
    {
        fputc('0', file);
    }
    else if (len == 1)
    {
        _fmpq_fprint(file, poly + 0, den);
    }
    else if (len == 2)
    {
        if (__fmpq_is_one(poly + 1, den))
        {
            flint_fprintf(file, "%s", x);
        }
        else if (__fmpq_is_pm1(poly + 1, den))
        {
            flint_fprintf(file, "-%s", x);
        }
        else
        {
            __fmpq_fprint(poly + 1, den);
            flint_fprintf(file, "*%s", x);
        }

        if (fmpz_sgn(poly + 0) > 0)
        {
            flint_fprintf(file, "+");
            __fmpq_fprint(poly + 0, den);
        }
        else if (fmpz_sgn(poly + 0) < 0)
        {
            __fmpq_fprint(poly + 0, den);
        }
    }
    else  /* len >= 3 */
    {
        slong i = len - 1;  /* i >= 2 */
        {
            if (__fmpq_is_one(poly + i, den))
               flint_fprintf(file, "%s^%wd", x, i);
            else if (__fmpq_is_pm1(poly + i, den))
               flint_fprintf(file, "-%s^%wd", x, i);
            else
            {
               __fmpq_fprint(poly + i, den);
               flint_fprintf(file, "*%s^%wd", x, i);
            }
            --i;
        }

        for (; i > 1; --i)
        {
            if (poly[i] == 0)
                continue;

            if (__fmpq_is_one(poly + i, den))
                flint_fprintf(file, "+%s^%wd", x, i);
            else if (__fmpq_is_pm1(poly + i, den))
                flint_fprintf(file, "-%s^%wd", x, i);
            else
            {
                if (fmpz_sgn(poly + i) > 0)
                {
                    fputc('+', file);
                }
                __fmpq_fprint(poly + i, den);
                flint_fprintf(file, "*%s^%wd", x, i);
            }
        }

        if (poly[1])
        {
            if (__fmpq_is_one(poly + 1, den))
            {
                fputc('+', file);
                fputs(x, file);
            }
            else if (__fmpq_is_pm1(poly + 1, den))
            {
                fputc('-', file);
                fputs(x, file);
            }
            else
            {
                if (fmpz_sgn(poly + 1) > 0)
                {
                    fputc('+', file);
                }
                __fmpq_fprint(poly + 1, den);
                fputc('*', file);
                fputs(x, file);
            }
        }
        if (*(poly))
        {
            if (fmpz_sgn(poly) > 0)
            {
                fputc('+', file);
            }
            __fmpq_fprint(poly + 0, den);
        }
    }

    fmpz_clear(n);
    fmpz_clear(d);
    fmpz_clear(g);

    return 1;
}

#undef __fmpq_fprint

int fmpq_poly_fprint(FILE * file, const fmpq_poly_t poly) { return _fmpq_poly_fprint(file, poly->coeffs, poly->den, poly->length); }
int fmpq_poly_fprint_pretty(FILE * file, const fmpq_poly_t poly, const char * var) { return _fmpq_poly_fprint_pretty(file, poly->coeffs, poly->den, poly->length, var); }
int _fmpq_poly_print(const fmpz * poly, const fmpz_t den, slong len) { return _fmpq_poly_fprint(stdout, poly, den, len); }
int fmpq_poly_print(const fmpq_poly_t poly) { return fmpq_poly_fprint(stdout, poly); }
int _fmpq_poly_print_pretty(const fmpz *poly, const fmpz_t den, slong len, const char * x) { return _fmpq_poly_fprint_pretty(stdout, poly, den, len, x); }
int fmpq_poly_print_pretty(const fmpq_poly_t poly, const char * var) { return fmpq_poly_fprint_pretty(stdout, poly, var); }

/* reading ********************************************************************/

int fmpq_poly_fread(FILE * file, fmpq_poly_t poly)
{
    int r;
    slong i, len;
    mpz_t t;
    mpq_t *a;

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

    a = flint_malloc(len * sizeof(mpq_t));
    for (i = 0; i < len; i++)
        mpq_init(a[i]);

    for (i = 0; (i < len) && r; i++)
        r = mpq_inp_str(a[i], file, 10);

    if (r > 0)
        __fmpq_poly_set_array_mpq(poly, (const mpq_t *) a, len);

    for (i = 0; i < len; i++)
        mpq_clear(a[i]);
    flint_free(a);

    return r;
}

int fmpq_poly_read(fmpq_poly_t poly) { return fmpq_poly_fread(stdin, poly); }

/* debugging ******************************************************************/

int fmpq_poly_debug(const fmpq_poly_t poly)
{
    slong i;

    flint_printf("{alloc: %wd, length: %wd, coeffs:", poly->alloc, poly->length);
    for (i = 0; i < poly->alloc; i++)
    {
        flint_printf(" ");
        fmpz_print(poly->coeffs + i);
    }
    flint_printf(", den: ");
    fmpz_print(poly->den);
    flint_printf("}");

    return 1;
}

