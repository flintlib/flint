/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#undef ulong
#define ulong ulongxx /* ensure vendor doesn't typedef ulong */
#include <stdio.h>
#undef ulong
#include "fmpq_mat.h"

void fmpq_mat_print(const fmpq_mat_t mat)
{
    slong i, j;

    printf("<" WORD_FMT "d x " WORD_FMT "d matrix over Q>\n", mat->r, mat->c);

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
