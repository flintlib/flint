/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include <stdlib.h>
#include "flint.h"
#include "ulong_extras.h"
#include "nmod_poly.h"

void
nmod_poly_randtest(nmod_poly_t poly, flint_rand_t state, slong len)
{
    nmod_poly_fit_length(poly, len);
    _nmod_vec_randtest(poly->coeffs, state, len, poly->mod);
    poly->length = len;
    _nmod_poly_normalise(poly);
}

void
nmod_poly_randtest_monic(nmod_poly_t poly, flint_rand_t state, slong len)
{
    nmod_poly_fit_length(poly, len);
    _nmod_vec_randtest(poly->coeffs, state, len - 1, poly->mod);
    poly->coeffs[len - 1] = 1;
    poly->length = len;
}

static void
nmod_poly_randtest_monic_sparse(nmod_poly_t poly, flint_rand_t state,
                                                      slong len, slong nonzero)
{
    slong i;

    nmod_poly_fit_length(poly, len);
    _nmod_vec_zero(poly->coeffs, len);
    poly->coeffs[0] = n_randtest(state) % poly->mod.n;
    for (i = 1; i < nonzero; i++)
       poly->coeffs[n_randint(state, len - 1) + 1] = n_randtest(state) % poly->mod.n;
    poly->coeffs[len - 1] = 1;
    _nmod_poly_set_length(poly, len);
}

void
nmod_poly_randtest_irreducible(nmod_poly_t poly, flint_rand_t state, slong len)
{
    do {
        nmod_poly_randtest(poly, state, len);
    } while (nmod_poly_is_zero(poly) || !(nmod_poly_is_irreducible(poly)));
}

void
nmod_poly_randtest_monic_irreducible(nmod_poly_t poly, flint_rand_t state, slong len)
{
    do {
        nmod_poly_randtest_monic(poly, state, len);
    } while (nmod_poly_is_zero(poly) || !(nmod_poly_is_irreducible(poly)));
}

static void
nmod_poly_randtest_monic_irreducible_sparse(nmod_poly_t poly,
                                              flint_rand_t state, slong len)
{
    slong i = 0;
    slong terms = 3;
    do {
        i++;
        terms += ((i % 4) == 0);
        if (terms >= len)
            terms = 3;
        nmod_poly_randtest_monic_sparse(poly, state, len, terms);
    } while (nmod_poly_is_zero(poly) ||
             !nmod_poly_is_irreducible(poly));
}

void
nmod_poly_randtest_trinomial(nmod_poly_t poly, flint_rand_t state, slong len)
{
    ulong k;
    nmod_poly_fit_length(poly, len);
    _nmod_vec_zero(poly->coeffs, len);
    poly->coeffs[0] = n_randtest(state) % poly->mod.n;
    poly->coeffs[len - 1] = 1;
    k = (n_randtest(state) % (len - 2)) + 1;
    poly->coeffs[k] = n_randtest(state) % poly->mod.n;
    _nmod_poly_set_length(poly, len);
}

void
nmod_poly_randtest_pentomial(nmod_poly_t poly, flint_rand_t state, slong len)
{
    nmod_poly_fit_length(poly, len);
    _nmod_vec_zero(poly->coeffs, len);
    poly->coeffs[0] = n_randtest(state) % poly->mod.n;
    poly->coeffs[1] = n_randtest(state) % poly->mod.n;
    poly->coeffs[2] = n_randtest(state) % poly->mod.n;
    poly->coeffs[3] = n_randtest(state) % poly->mod.n;
    poly->coeffs[len - 1] = 1;
    _nmod_poly_set_length(poly, len);
}

int
nmod_poly_randtest_trinomial_irreducible(nmod_poly_t poly, flint_rand_t state,
                                         slong len, slong max_attempts)
{
    slong i = 0;

    while (max_attempts == 0 || i < max_attempts)
    {
        nmod_poly_randtest_trinomial(poly, state, len);
        if (!nmod_poly_is_zero(poly) && nmod_poly_is_irreducible(poly))
        {
            return 1;
        }
        i++;
        
    }
    return 0;
}

int
nmod_poly_randtest_pentomial_irreducible(nmod_poly_t poly, flint_rand_t state,
                                         slong len, slong max_attempts)
{
    slong i = 0;

    while (max_attempts == 0 || i < max_attempts)
    {
        nmod_poly_randtest_pentomial(poly, state, len);
        if (!nmod_poly_is_zero(poly) && nmod_poly_is_irreducible(poly))
        {
            return 1;
        }
        i++;
        
    }
    return 0;
}

void
nmod_poly_randtest_sparse_irreducible(nmod_poly_t poly, flint_rand_t state, slong len)
{
    if (len < 3)
    {
        nmod_poly_randtest_monic_irreducible(poly, state, len);
        return;
    }

    /* Try trinomials */
    if (nmod_poly_randtest_trinomial_irreducible(poly, state, len, 2*len))
        return;

    if (len < 5)
    {
        nmod_poly_randtest_monic_irreducible(poly, state, len);
        return;
    }

    /* Try pentomials */
    if (nmod_poly_randtest_pentomial_irreducible(poly, state, len, 2*len))
        return;

    /* Random monic sparse */
    nmod_poly_randtest_monic_irreducible_sparse(poly, state, len);
}
