/*
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2009 William Hart
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"
#include "fmpz_mod_poly_factor.h"

void fmpz_mod_poly_randtest(fmpz_mod_poly_t f, flint_rand_t state, slong len,
                                                      const fmpz_mod_ctx_t ctx)
{
    slong i;

    fmpz_mod_poly_fit_length(f, len, ctx);

    for (i = 0; i < len; i++)
        fmpz_randm(f->coeffs + i, state, fmpz_mod_ctx_modulus(ctx));

    _fmpz_mod_poly_set_length(f, len);
    _fmpz_mod_poly_normalise(f);
}

void fmpz_mod_poly_randtest_monic(fmpz_mod_poly_t f, flint_rand_t state,
                                           slong len, const fmpz_mod_ctx_t ctx)
{
    slong i;

    FLINT_ASSERT(len > 0);

    fmpz_mod_poly_fit_length(f, len, ctx);

    for (i = 0; i < len - 1; i++)
        fmpz_randm(f->coeffs + i, state, fmpz_mod_ctx_modulus(ctx));

    fmpz_one(f->coeffs + len - 1);

    _fmpz_mod_poly_set_length(f, len);
}

static void
fmpz_mod_poly_randtest_monic_sparse(fmpz_mod_poly_t poly, flint_rand_t state,
                            slong len, slong nonzero, const fmpz_mod_ctx_t ctx)
{
    slong i;

    fmpz_mod_poly_fit_length(poly, len, ctx);
    _fmpz_vec_zero(poly->coeffs, len);
    fmpz_randm(poly->coeffs + 0, state, fmpz_mod_ctx_modulus(ctx));
    for (i = 1; i < nonzero; i++)
       fmpz_randm(poly->coeffs + n_randint(state, len - 1) + 1,
                                             state, fmpz_mod_ctx_modulus(ctx));
    fmpz_set_ui(poly->coeffs + len - 1, 1);
    _fmpz_mod_poly_set_length(poly, len);
}

void fmpz_mod_poly_randtest_irreducible(fmpz_mod_poly_t f, flint_rand_t state,
                                           slong len, const fmpz_mod_ctx_t ctx)
{
    if (len == 0)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_mod_poly_randtest_irreducible). len == 0.\n");
    }

    do {
        fmpz_mod_poly_randtest(f, state, len, ctx);
    } while (fmpz_mod_poly_is_zero(f, ctx) ||
             !fmpz_mod_poly_is_irreducible(f, ctx));
}

void fmpz_mod_poly_randtest_monic_irreducible(fmpz_mod_poly_t f,
                      flint_rand_t state, slong len, const fmpz_mod_ctx_t ctx)
{
    if (len == 0)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_mod_poly_randtest_monic_irreducible). len == 0.\n");
    }

    do {
        fmpz_mod_poly_randtest_monic(f, state, len, ctx);
    } while (fmpz_mod_poly_is_zero(f, ctx) ||
             !fmpz_mod_poly_is_irreducible(f, ctx));
}

void fmpz_mod_poly_randtest_not_zero(fmpz_mod_poly_t f,
                       flint_rand_t state, slong len, const fmpz_mod_ctx_t ctx)
{
    if (len == 0)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_mod_poly_randtest_not_zero). len == 0.\n");
    }

    do {
        fmpz_mod_poly_randtest(f, state, len, ctx);
    } while (fmpz_mod_poly_is_zero(f, ctx));
}

static void
fmpz_mod_poly_randtest_monic_irreducible_sparse(fmpz_mod_poly_t poly,
                       flint_rand_t state, slong len, const fmpz_mod_ctx_t ctx)
{
    slong i = 0;
    slong terms = 3;
    do {
        i++;
        terms += ((i % 4) == 0);
        if (terms >= len)
            terms = 3;
        fmpz_mod_poly_randtest_monic_sparse(poly, state, len, terms, ctx);
    } while (fmpz_mod_poly_is_zero(poly, ctx) ||
             !fmpz_mod_poly_is_irreducible(poly, ctx));
}

void fmpz_mod_poly_randtest_trinomial(fmpz_mod_poly_t poly,
                       flint_rand_t state, slong len, const fmpz_mod_ctx_t ctx)
{
    ulong k;
    fmpz_mod_poly_fit_length(poly, len, ctx);
    _fmpz_vec_zero(poly->coeffs, len);
    fmpz_randm(poly->coeffs, state, fmpz_mod_ctx_modulus(ctx));
    k = (n_randtest(state) % (len - 2)) + 1;
    fmpz_randm(poly->coeffs + k, state, fmpz_mod_ctx_modulus(ctx));
    fmpz_one(poly->coeffs + len - 1);
    _fmpz_mod_poly_set_length(poly, len);
}

void fmpz_mod_poly_randtest_pentomial(fmpz_mod_poly_t poly,
                       flint_rand_t state, slong len, const fmpz_mod_ctx_t ctx)
{
    fmpz_mod_poly_fit_length(poly, len, ctx);
    _fmpz_vec_zero(poly->coeffs, len);
    fmpz_randm(poly->coeffs, state, fmpz_mod_ctx_modulus(ctx));
    fmpz_randm(poly->coeffs + 1, state, fmpz_mod_ctx_modulus(ctx));
    fmpz_randm(poly->coeffs + 2, state, fmpz_mod_ctx_modulus(ctx));
    fmpz_randm(poly->coeffs + 3, state, fmpz_mod_ctx_modulus(ctx));
    fmpz_one(poly->coeffs + len - 1);
    _fmpz_mod_poly_set_length(poly, len);
}

int fmpz_mod_poly_randtest_trinomial_irreducible(fmpz_mod_poly_t poly,
                           flint_rand_t state, slong len, slong max_attempts,
                                                      const fmpz_mod_ctx_t ctx)
{
    slong i = 0;

    while (max_attempts == 0 || i < max_attempts)
    {
        fmpz_mod_poly_randtest_trinomial(poly, state, len, ctx);
        if (!fmpz_mod_poly_is_zero(poly, ctx) &&
            fmpz_mod_poly_is_irreducible(poly, ctx))
        {
            return 1;
        }
        i++;
    }
    return 0;
}

int fmpz_mod_poly_randtest_pentomial_irreducible(fmpz_mod_poly_t poly,
                           flint_rand_t state, slong len, slong max_attempts,
                                                      const fmpz_mod_ctx_t ctx)
{
    slong i = 0;

    while (max_attempts == 0 || i < max_attempts)
    {
        fmpz_mod_poly_randtest_pentomial(poly, state, len, ctx);
        if (!fmpz_mod_poly_is_zero(poly, ctx) &&
            fmpz_mod_poly_is_irreducible(poly, ctx))
        {
            return 1;
        }
        i++;
    }
    return 0;
}

void fmpz_mod_poly_randtest_sparse_irreducible(fmpz_mod_poly_t poly,
                       flint_rand_t state, slong len, const fmpz_mod_ctx_t ctx)
{
    if (len < 3)
    {
        fmpz_mod_poly_randtest_monic_irreducible(poly, state, len, ctx);
        return;
    }

    /* Try trinomials */
    if (fmpz_mod_poly_randtest_trinomial_irreducible(poly, state, len, 2*len, ctx))
        return;

    if (len < 5)
    {
        fmpz_mod_poly_randtest_monic_irreducible(poly, state, len, ctx);
        return;
    }

    /* Try pentomials */
    if (fmpz_mod_poly_randtest_pentomial_irreducible(poly, state, len, 2*len, ctx))
        return;

    /* Random monic sparse */
    fmpz_mod_poly_randtest_monic_irreducible_sparse(poly, state, len, ctx);
}
