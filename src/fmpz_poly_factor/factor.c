/*
    Copyright (C) 2011 Andy Novocin
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_poly_factor.h"

static void fmpz_poly_factor_deflation(fmpz_poly_factor_t fac,
                                       const fmpz_poly_t G, int deflation);

static void _fmpz_poly_factor_inflation(fmpz_poly_factor_t fac,
                                        const fmpz_poly_t T, ulong d,
                                        slong exp);

static ulong
_smallest_prime_factor(ulong d)
{
    ulong t;
    if (d % 2 == 0)
        return 2;
    for (t = 3; t * t <= d; t += 2)
        if (d % t == 0)
            return t;
    return d;
}

/*
    Append the factorisation of T(x^d), with outer multiplicity exp, to fac,
    where T is irreducible over Q, primitive with positive leading
    coefficient.  Chains over the prime factors p of d via
    T(x^d) = U(x^(d/p)) with U = T(x^p): each prime step is either certified
    irreducible by _fmpz_poly_factor_inflation_is_irreducible_capelli
    (in which case no factorisation
    is performed at all) or handled by the generic machinery at degree
    p*deg(T), which is never larger, and usually much smaller, than the
    single call at degree d*deg(T) it replaces.
*/
static void
_fmpz_poly_factor_inflation(fmpz_poly_factor_t fac, const fmpz_poly_t T,
                            ulong d, slong exp)
{
    ulong p;
    fmpz_poly_t U;

    if (d == 1)
    {
        fmpz_poly_factor_insert(fac, T, exp);
        return;
    }

    /* T(0) != 0 is guaranteed: fmpz_poly_factor_deflation strips powers
       of x before deflating, and factors of U = T(x^p) with U(0) != 0
       again have nonzero constant term, so the invariant persists through
       the recursion below. */

    p = _smallest_prime_factor(d);

    fmpz_poly_init(U);
    fmpz_poly_inflate(U, T, p);

    if (_fmpz_poly_factor_inflation_is_irreducible_capelli(T, p))
    {
        /* U = T(x^p) is certified irreducible */
        _fmpz_poly_factor_inflation(fac, U, d / p, exp);
    }
    else
    {
        slong i;
        fmpz_poly_factor_t ufac;

        fmpz_poly_factor_init(ufac);
        fmpz_poly_factor_deflation(ufac, U, 0);

        for (i = 0; i < ufac->num; i++)
            _fmpz_poly_factor_inflation(fac, ufac->p + i, d / p,
                                        exp * ufac->exp[i]);

        fmpz_poly_factor_clear(ufac);
    }

    fmpz_poly_clear(U);
}

static void fmpz_poly_factor_deflation(fmpz_poly_factor_t fac, const fmpz_poly_t G, int deflation)
{
    const slong lenG = G->length;
    fmpz_poly_t g;

    fac->num = 0;

    if (lenG <= 1)
    {
        if (lenG < 1)
            fmpz_zero(&fac->c);
        else
            fmpz_set(&fac->c, G->coeffs + 0);
        return;
    }

    fmpz_poly_init(g);

    if (lenG < 5)
    {
        fmpz_poly_content(&fac->c, G);
        if (fmpz_sgn(fmpz_poly_lead(G)) < 0)
            fmpz_neg(&fac->c, &fac->c);
        fmpz_poly_scalar_divexact_fmpz(g, G, &fac->c);

        if (lenG < 3)
            fmpz_poly_factor_insert(fac, g, 1);
        else if (lenG == 3)
            _fmpz_poly_factor_quadratic(fac, g, 1);
        else
            _fmpz_poly_factor_cubic(fac, g, 1);
    }
    else
    {
        slong i, j, k, d;
        fmpz_poly_factor_t sq_fr_fac;

        /* Does a presearch for a factor of form x^k */
        for (k = 0; fmpz_is_zero(G->coeffs + k); k++) ;

        if (k != 0)
        {
            fmpz_poly_t t;

            fmpz_poly_init(t);
            fmpz_poly_set_coeff_ui(t, 1, 1);
            fmpz_poly_factor_insert(fac, t, k);
            fmpz_poly_clear(t);
        }

        fmpz_poly_shift_right(g, G, k);

        if (deflation && (d = fmpz_poly_deflation(G)) > 1)
        {
            fmpz_poly_factor_t gfac;
            fmpz_poly_factor_init(gfac);

            fmpz_poly_deflate(g, g, d);
            fmpz_poly_factor(gfac, g);
            fmpz_set(&fac->c, &gfac->c);

            /* Each gfac->p + i is irreducible; factor its inflation by
               chained Capelli certification instead of refactoring
               (gfac->p + i)(x^d) from scratch. */
            for (i = 0; i < gfac->num; i++)
                _fmpz_poly_factor_inflation(fac, gfac->p + i, (ulong) d,
                                            gfac->exp[i]);

            fmpz_poly_factor_clear(gfac);
        }
        else
        {
            /* Could make other tests for x-1 or simple things
               maybe take advantage of the composition algorithm */
            fmpz_poly_factor_init(sq_fr_fac);
            fmpz_poly_factor_squarefree(sq_fr_fac, g);

            fmpz_set(&fac->c, &sq_fr_fac->c);

            /* Factor each square-free part */
            for (j = 0; j < sq_fr_fac->num; j++)
            {
                _fmpz_poly_factor_zassenhaus(fac, sq_fr_fac->exp[j],
                                                           sq_fr_fac->p + j, 8, 1);
            }

            fmpz_poly_factor_clear(sq_fr_fac);
        }
    }
    fmpz_poly_clear(g);
}

void fmpz_poly_factor(fmpz_poly_factor_t fac, const fmpz_poly_t G)
{
    fmpz_poly_factor_deflation(fac, G, 1);
}
