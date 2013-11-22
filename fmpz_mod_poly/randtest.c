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
    Copyright (C) 2013 Mike Hansen

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mod_poly.h"

void
fmpz_mod_poly_randtest(fmpz_mod_poly_t f, flint_rand_t state, slong len)
{
    slong i;

    fmpz_mod_poly_fit_length(f, len);

    for (i = 0; i < len; i++)
        fmpz_randm(f->coeffs + i, state, &(f->p));

    _fmpz_mod_poly_set_length(f, len);
    _fmpz_mod_poly_normalise(f);
}

void
fmpz_mod_poly_randtest_monic(fmpz_mod_poly_t f, flint_rand_t state, slong len)
{
    slong i;

    fmpz_mod_poly_fit_length(f, len);

    for (i = 0; i < len - 1; i++)
        fmpz_randm(f->coeffs + i, state, &(f->p));

    fmpz_one(f->coeffs + len - 1);

    _fmpz_mod_poly_set_length(f, len);
}

void
fmpz_mod_poly_randtest_irreducible(fmpz_mod_poly_t f,
                                   flint_rand_t state, slong len)
{
    if (len == 0)
    {
        flint_printf("Exception (fmpz_mod_poly_randtest_irreducible). len == 0.\n");
        abort();
    }

    fmpz_mod_poly_randtest(f, state, len);
    while (fmpz_mod_poly_is_zero(f) || !fmpz_mod_poly_is_irreducible(f))
        fmpz_mod_poly_randtest(f, state, len);
}

void
fmpz_mod_poly_randtest_monic_irreducible(fmpz_mod_poly_t f,
                                         flint_rand_t state, slong len)
{
    if (len == 0)
    {
        flint_printf("Exception (fmpz_mod_poly_randtest_monic_irreducible). len == 0.\n");
        abort();
    }

    fmpz_mod_poly_randtest_monic(f, state, len);
    while (fmpz_mod_poly_is_zero(f) || !fmpz_mod_poly_is_irreducible(f))
        fmpz_mod_poly_randtest_monic(f, state, len);
}

void
fmpz_mod_poly_randtest_not_zero(fmpz_mod_poly_t f,
                                flint_rand_t state, slong len)
{
    if (len == 0)
    {
        flint_printf("Exception (fmpz_mod_poly_randtest_not_zero). len == 0.\n");
        abort();
    }

    fmpz_mod_poly_randtest(f, state, len);
    while (fmpz_mod_poly_is_zero(f))
        fmpz_mod_poly_randtest(f, state, len);
}

void
fmpz_mod_poly_randtest_trinomial(fmpz_mod_poly_t poly, flint_rand_t state, slong len)
{
    ulong k;
    fmpz_mod_poly_fit_length(poly, len);
    _fmpz_vec_zero(poly->coeffs, len);
    fmpz_randm(poly->coeffs, state, &(poly->p));
    k = (n_randtest(state) % (len - 2)) + 1;
    fmpz_randm(poly->coeffs + k, state, &(poly->p));
    fmpz_one(poly->coeffs + len - 1);
    _fmpz_mod_poly_set_length(poly, len);
}

void
fmpz_mod_poly_randtest_pentomial(fmpz_mod_poly_t poly, flint_rand_t state, slong len)
{
    fmpz_mod_poly_fit_length(poly, len);
    _fmpz_vec_zero(poly->coeffs, len);
    fmpz_randm(poly->coeffs, state, &(poly->p));
    fmpz_randm(poly->coeffs + 1, state, &(poly->p));
    fmpz_randm(poly->coeffs + 2, state, &(poly->p));
    fmpz_randm(poly->coeffs + 3, state, &(poly->p));
    fmpz_one(poly->coeffs + len - 1);
    _fmpz_mod_poly_set_length(poly, len);
}

int
fmpz_mod_poly_randtest_trinomial_irreducible(fmpz_mod_poly_t poly, flint_rand_t state,
                                             slong len, slong max_attempts)
{
    slong i = 0;

    while (max_attempts == 0 || i < max_attempts)
    {
        fmpz_mod_poly_randtest_trinomial(poly, state, len);
        if (!fmpz_mod_poly_is_zero(poly) && fmpz_mod_poly_is_irreducible(poly))
        {
            return 1;
        }
        i++;
        
    }
    return 0;
}

int
fmpz_mod_poly_randtest_pentomial_irreducible(fmpz_mod_poly_t poly, flint_rand_t state,
                                             slong len, slong max_attempts)
{
    slong i = 0;

    while (max_attempts == 0 || i < max_attempts)
    {
        fmpz_mod_poly_randtest_pentomial(poly, state, len);
        if (!fmpz_mod_poly_is_zero(poly) && fmpz_mod_poly_is_irreducible(poly))
        {
            return 1;
        }
        i++;
        
    }
    return 0;
}

void
fmpz_mod_poly_randtest_sparse_irreducible(fmpz_mod_poly_t poly, flint_rand_t state, slong len)
{
    if (len < 3)
    {
        fmpz_mod_poly_randtest_monic_irreducible(poly, state, len);
        return;
    }

    /* Try trinomials */
    if (fmpz_mod_poly_randtest_trinomial_irreducible(poly, state, len, 2*len))
        return;

    if (len < 5)
    {
        fmpz_mod_poly_randtest_monic_irreducible(poly, state, len);
        return;
    }

    /* Try pentomials */
    if (fmpz_mod_poly_randtest_pentomial_irreducible(poly, state, len, 2*len))
        return;

    /* Give up */
    fmpz_mod_poly_randtest_monic_irreducible(poly, state, len);
}
