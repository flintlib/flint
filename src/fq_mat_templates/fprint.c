/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

/*
    The macros xxx_putc, xxx_printf, and xxx_T_print are provided as
    wrappers to handle return values and error conditions.  While this
    is not exactly pretty, it improves the readability of the
    functions TEMPLATE(T, mat_fprint) and TEMPLATE(T,
    mat_fprint_pretty).  Moreover, if we later want to improve the
    handling of returns values, e.g.  to return the number of
    characters printed, this will be easier.

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
    z = flint_fprintf(file, WORD_FMT "d " WORD_FMT "d  ", r, c);  \
    if (z <= 0)                            \
        return z;                          \
} while (0)

#define xxx_T_print(f, ctx)                     \
do {                                            \
    z = TEMPLATE(T, fprint)(file, (f), (ctx));  \
    if (z <= 0)                                 \
        return z;                               \
} while(0)

#define xxx_T_print_pretty(f, ctx)                     \
do {                                                   \
    z = TEMPLATE(T, fprint_pretty)(file, (f), (ctx));  \
    if (z <= 0)                                        \
        return z;                                      \
} while(0)

int TEMPLATE(T, mat_fprint) (FILE * file, const TEMPLATE(T, mat_t) mat,
                             const TEMPLATE(T, ctx_t) ctx)
{
    int z;
    slong i, j;
    slong r = mat->r;
    slong c = mat->c;

    xxx_printf();
    for (i = 0; (i < r); i++)
    {
        for (j = 0; j < c; j++)
        {
            xxx_T_print(mat->rows[i] + j, ctx);
            if (j != c - 1)
                xxx_putc(' ');
        }
        if (i != r - 1)
            xxx_putc(' ');
    }

    return z;
}

int TEMPLATE(T, mat_fprint_pretty) (FILE * file, const TEMPLATE(T, mat_t) mat,
                                    const TEMPLATE(T, ctx_t) ctx)
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
            xxx_T_print_pretty(mat->rows[i] + j, ctx);
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
#undef xxx_T_print

#endif
