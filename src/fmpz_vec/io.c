/*
    Copyright (C) 2008, 2009, 2010 William Hart
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "gmpcompat.h"
#include "fmpz.h"
#include "fmpz_vec.h"

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

int _fmpz_vec_fprint(FILE * file, const fmpz * vec, slong len)
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
                r = fmpz_fprint(file, vec + i);
        }
    }

    return r;
}

int _fmpz_vec_print(const fmpz * vec, slong len) { return _fmpz_vec_fprint(stdout, vec, len); }

/* reading ********************************************************************/

int _fmpz_vec_fread(FILE * file, fmpz ** vec, slong * len)
{
    int alloc, r;
    slong i;
    mpz_t t;

    alloc = (*vec == NULL);

    mpz_init(t);
    r = mpz_inp_str(t, file, 10);
    if (r == 0)
    {
        if (alloc)
            *len = 0;
        mpz_clear(t);
        return 0;
    }
    if (!mpz_fits_slong_p(t))
    {
        flint_throw(FLINT_ERROR, "Exception (_fmpz_vec_fread). Length does not fit into a slong.\n");
    }
    if (alloc)
    {
        *len = flint_mpz_get_si(t);
        *vec = _fmpz_vec_init(*len);
    }
    else
    {
        if (*len != flint_mpz_get_si(t))
        {
            mpz_clear(t);
            return 0;
        }
    }
    mpz_clear(t);

    for (i = 0; i < *len; i++)
    {
        r = fmpz_fread(file, (*vec) + i);
        if (r <= 0)
        {
            if (alloc)
            {
                _fmpz_vec_clear(*vec, *len);
                *vec = NULL;
                *len = 0;
            }
            return r;
        }
    }

    return 1;
}

int _fmpz_vec_read(fmpz ** vec, slong * len) { return _fmpz_vec_fread(stdin, vec, len); }
