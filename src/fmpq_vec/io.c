/*
    Copyright (C) 2016 Vincent Delecroix

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "fmpq.h"
#include "fmpq_vec.h"

/* printing *******************************************************************/

int _fmpq_vec_fprint(FILE * file, const fmpq * vec, slong len)
{
    int r;
    slong i;

    r = flint_fprintf(file, "%wd", len);
    if ((len > 0) && (r > 0))
    {
        r = fputc(' ', file);
        for (i = 0; (i < len) && (r > 0); i++)
        {
            r = fputc(' ', file);
            if (r > 0)
                r = fmpq_fprint(file, vec + i);
        }
    }

    return r;
}

int _fmpq_vec_print(const fmpq * vec, slong len) { return _fmpq_vec_fprint(stdout, vec, len); }
