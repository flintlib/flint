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

    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#include <stdio.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"

int fmpq_poly_fread(FILE * file, fmpq_poly_t poly)
{
    int r;
    long i, len;
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
        printf("Exception (fmpz_poly_fread). Length does not fit into a long.\n");
        abort();
    }
    len = mpz_get_si(t);
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

