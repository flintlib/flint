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

    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2009 William Hart

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mod_poly.h"
#include "fmpz_mod_poly_factor.h"

void
fmpz_mod_poly_randtest(fmpz_mod_poly_t f, flint_rand_t state, long len)
{
    long i;

    fmpz_mod_poly_fit_length(f, len);

    for (i = 0; i < len; i++)
        fmpz_randm(f->coeffs + i, state, &(f->p));

    _fmpz_mod_poly_set_length(f, len);
    _fmpz_mod_poly_normalise(f);
}

void
fmpz_mod_poly_randtest_irreducible(fmpz_mod_poly_t f,
                                   flint_rand_t state, long len)
{
    if (len == 0)
    {
        printf("Exception (fmpz_mod_poly_randtest_irreducible). len == 0.\n");
        abort();
    }

    fmpz_mod_poly_randtest(f, state, len);
    while (fmpz_mod_poly_is_zero(f) || !fmpz_mod_poly_is_irreducible(f))
        fmpz_mod_poly_randtest(f, state, len);
}

void
fmpz_mod_poly_randtest_not_zero(fmpz_mod_poly_t f,
                                flint_rand_t state, long len)
{
    if (len == 0)
    {
        printf("Exception (fmpz_mod_poly_randtest_not_zero). len == 0.\n");
        abort();
    }

    fmpz_mod_poly_randtest(f, state, len);
    while (fmpz_mod_poly_is_zero(f))
        fmpz_mod_poly_randtest(f, state, len);
}

