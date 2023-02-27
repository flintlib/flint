/*
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"

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
        flint_printf("Exception (fmpz_poly_fread). Length does not fit into a slong.\n");
        flint_abort();
    }
    len = flint_mpz_get_si(t);
    mpz_clear(t);

    a = flint_malloc(len * sizeof(mpq_t));
    for (i = 0; i < len; i++)
        mpq_init(a[i]);

    for (i = 0; (i < len) && r; i++)
        r = mpq_inp_str(a[i], file, 10);

    if (r > 0)
        fmpq_poly_set_array_mpq(poly, (const mpq_t *) a, len);

    for (i = 0; i < len; i++)
        mpq_clear(a[i]);
    flint_free(a);

    return r;
}

