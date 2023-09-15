/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "fmpq_mat.h"

/* printing *******************************************************************/

void fmpq_mat_print(const fmpq_mat_t mat)
{
    slong i, j;

    flint_printf("<%wd x %wd matrix over Q>\n", mat->r, mat->c);

    for (i = 0; i < mat->r; i++)
    {
        flint_printf("[");
        for (j = 0; j < mat->c; j++)
        {
            fmpq_print(fmpq_mat_entry(mat, i, j));
            if (j + 1 < mat->c)
                flint_printf(", ");
        }
        flint_printf("]\n");
    }
    flint_printf("\n");
}


/* This is copied from flint/src/fmpz_mat/io.c */

#ifdef FLINT_HAVE_FILE

#define xxx_putc(c)        \
do {                       \
    z = fputc((c), file);  \
    if (z <= 0)            \
        return z;          \
} while (0)

#define xxx_fmpq_print(f)        \
do {                             \
    z = fmpq_fprint(file, (f));  \
    if (z <= 0)                  \
        return z;                \
} while(0)

int fmpq_mat_fprint_pretty(FILE * file, const fmpq_mat_t mat)
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
            xxx_fmpq_print(mat->rows[i] + j);
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
#undef xxx_fmpq_print

#endif
