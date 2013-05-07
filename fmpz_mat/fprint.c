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

    Copyright (C) 2010 William Hart

******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mat.h"

/*
    The macros xxx_putc, xxx_printf, and xxx_fmpz_print are provided 
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

#define xxx_printf()                       \
do {                                       \
    z = fprintf(file, "%li %li  ", r, c);  \
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
    long i, j;
    long r = mat->r;
    long c = mat->c;

    xxx_printf();
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
    long i, j;
    long r = mat->r;
    long c = mat->c;

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
#undef xxx_printf
#undef xxx_fmpz_print

