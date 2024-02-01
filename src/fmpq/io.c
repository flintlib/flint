/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "fmpq.h"

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

int _fmpq_fprint(FILE * file, const fmpz_t num, const fmpz_t den)
{
    if (fmpz_is_one(den))
    {
        return fmpz_fprint(file, num);
    }
    else
    {
        int r;
        r = fmpz_fprint(file, num);
        if (r > 0)
        {
            r = fputc('/', file);
            if (r > 0)
                r = fmpz_fprint(file, den);
        }
        return r;
    }
}

int fmpq_fprint(FILE * file, const fmpq_t x) { return _fmpq_fprint(file, &x->num, &x->den); }
int _fmpq_print(const fmpz_t num, const fmpz_t den) { return _fmpq_fprint(stdout, num, den); }
int fmpq_print(const fmpq_t x) { return fmpq_fprint(stdout, x); }
