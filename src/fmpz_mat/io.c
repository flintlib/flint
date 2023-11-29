/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2011  Andy Novocin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "gmpcompat.h"
#include "fmpz.h"
#include "fmpz_mat.h"

/* printing *******************************************************************/

/*
    The macros xxx_putc, xxx_flint_printf, and xxx_fmpz_print are provided
    as wrappers to handle return values and error conditions.  While
    this is not exactly pretty, it improves the readability of the
    functions fmpz_mat_fprint and fmpz_mat_fprint_pretty.  Moreover,
    if we later want to improve the handling of returns values, e.g.
    to return the number of characters printed, this will be easier.

    The macros are undef'd at the end of the file.
 */

#define xxx_putc(c)        \
do {                       \
    z = fputc((c), file);  \
    if (z <= 0)            \
        return z;          \
} while (0)

#define xxx_flint_printf()                       \
do {                                       \
    z = flint_fprintf(file, "%wd %wd  ", r, c);  \
    if (z <= 0)                            \
        return z;                          \
} while (0)

#define xxx_fmpz_print(f)        \
do {                             \
    z = fmpz_fprint(file, (f));  \
    if (z <= 0)                  \
        return z;                \
} while(0)

int fmpz_mat_fprint(FILE * file, const fmpz_mat_t mat)
{
    int z;
    slong i, j;
    slong r = mat->r;
    slong c = mat->c;

    xxx_flint_printf();
    for (i = 0; (i < r); i++)
    {
        for (j = 0; j < c; j++)
        {
            xxx_fmpz_print(mat->rows[i] + j);
            if (j != c - 1)
                xxx_putc(' ');
        }
        if (i != r - 1)
            xxx_putc(' ');
    }

    return z;
}

int fmpz_mat_fprint_pretty(FILE * file, const fmpz_mat_t mat)
{
    int z;
    slong i, j;
    slong r = mat->r;
    slong c = mat->c;

    xxx_putc('[');
    for (i = 0; i < r; i++)
    {
        xxx_putc('[');
        for (j = 0; j < c; j++)
        {
            xxx_fmpz_print(mat->rows[i] + j);
            if (j != c - 1)
                xxx_putc(' ');
        }
        xxx_putc(']');
        xxx_putc('\n');
    }
    xxx_putc(']');

    return z;
}

#undef xxx_putc
#undef xxx_flint_printf
#undef xxx_fmpz_print

int fmpz_mat_print(const fmpz_mat_t mat) { return fmpz_mat_fprint(stdout, mat); }
int fmpz_mat_print_pretty(const fmpz_mat_t mat) { return fmpz_mat_fprint_pretty(stdout, mat); }

/* reading ********************************************************************/

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
        flint_throw(FLINT_ERROR, "(fmpz_mat_fread): "
               "Number of rows does not fit into a slong.\n");
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
        flint_throw(FLINT_ERROR, "(fmpz_mat_fread): "
               "Number of columns does not fit into a slong.\n");
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
        flint_throw(FLINT_ERROR, "(fmpz_mat_fread): "
               "Dimensions are non-zero and do not match input dimensions.\n");
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

int fmpz_mat_read(fmpz_mat_t mat) { return fmpz_mat_fread(stdin, mat); }
