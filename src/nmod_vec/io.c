/*
    Copyright (C) 2023 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "ulong_extras.h"
#include "nmod_vec.h"

int _nmod_vec_fprint_pretty(FILE * file, nn_srcptr vec, slong len, nmod_t mod)
{
    slong j;
    int z, width;
    char fmt[FLINT_BITS + 5];

    z = flint_fprintf(file, "<length-%wd integer vector mod %wu>\n", len, mod.n);
    if (z <= 0)
        return z;

    if (!len)
        return z;

    width = n_sizeinbase(mod.n, 10);

    z = flint_sprintf(fmt, "%%%dwu", width);
    if (z <= 0)
        return z;

    z = flint_printf("[");
    if (z <= 0)
        return z;

    for (j = 0; j < len; j++)
    {
        z = flint_printf(fmt, vec[j]);
        if (z <= 0)
            return z;

        if (j + 1 < len)
        {
            z = flint_printf(" ");
            if (z <= 0)
                return z;
        }
    }

    z = flint_printf("]\n");

    return z;
}

void _nmod_vec_print_pretty(nn_srcptr vec, slong len, nmod_t mod)
{
    _nmod_vec_fprint_pretty(stdout, vec, len, mod);
}

int _nmod_vec_print(nn_srcptr vec, slong len, nmod_t mod)
{
    return _nmod_vec_fprint_pretty(stdout, vec, len, mod);
}

int _nmod_vec_fprint(FILE * f, nn_srcptr vec, slong len, nmod_t mod)
{
    return _nmod_vec_fprint_pretty(f, vec, len, mod);
}
