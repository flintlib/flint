/*
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "fmpz_poly.h"
#include "fmpz_poly_factor.h"
#include "arb_fmpz_poly.h"
#include "acb_poly.h"
#include "qqbar.h"

/*
We convert an annihilating polynomial over QQbar to an annihilating polynomial
over Z by multiplying out all combinations of conjugates of the coefficients.

TODO:

* Consider actually allowing squareful polynomials. In general, an initial
  squarefree factorization should be more efficient, but we *could* support
  squareful polynomials also here with a few adjustments.

* Consider using the quadratic formula for quadratics.

* Compute tighter denominator bounds and/or attempt to rescale things to get
  smaller coefficients.

* Exploit when some coefficients lie in a small common number field to attempt
  to reduce the degree.

* Is it ever faster to verify roots with an exact computation instead of
  numerics?

* Guess some low-degree roots with LLL and remove them before running the full
  factoring algorithm.

* Compute tight precision estimates.

*/

#define VERBOSE 0

/* Binary splitting product over all combinations of conjugates. */
static void
bsprod(acb_poly_t res, slong a, slong b, slong deg, qqbar_srcptr coeffs, acb_ptr * conj_roots, slong prec)
{
    if (b - a == 1)
    {
        slong k, j;
        acb_poly_fit_length(res, deg + 1);
        _acb_poly_set_length(res, deg + 1);
        k = a;
        for (j = 0; j <= deg; j++)
        {
            acb_set(res->coeffs + j, conj_roots[j] + (k % qqbar_degree(coeffs + j)));
            k /= qqbar_degree(coeffs + j);
        }
    }
    else
    {
        acb_poly_t f, g;

        acb_poly_init(f);
        acb_poly_init(g);

        bsprod(f, a, a + (b - a) / 2, deg, coeffs, conj_roots, prec);
        bsprod(g, a + (b - a) / 2, b, deg, coeffs, conj_roots, prec);
        acb_poly_mul(res, f, g, prec);

        acb_poly_clear(f);
        acb_poly_clear(g);
    }
}

int
_qqbar_roots_poly_squarefree(qqbar_ptr roots, qqbar_srcptr coeffs, slong len, slong deg_limit, slong bits_limit)
{
    slong deg, deg_bound, conjugations, i;
    slong prec, initial_prec;
    acb_ptr * conj_roots;
    acb_poly_t cpoly, ct;
    fmpz_poly_t rpoly;
    fmpz_t den, den_bound;
    int success = 1;

    deg = len - 1;

    if (deg_limit < 0)
        deg_limit = WORD_MAX;

    if (bits_limit < 0)
        bits_limit = WORD_MAX / 8;

    if (deg <= 0)
        return 1;

    FLINT_ASSERT(!qqbar_is_zero(coeffs + deg));

    if (deg == 1)
    {
        /* todo: check limits */
        qqbar_div(roots, coeffs, coeffs + 1);
        qqbar_neg(roots, roots);
        return 1;
    }

    conjugations = 1;
    for (i = 0; i <= deg; i++)
    {
        ulong hi, lo;
        umul_ppmm(hi, lo, conjugations, qqbar_degree(coeffs + i));
        if (hi != 0 || lo > (ulong) deg_limit)
            return 0;
        conjugations = lo;
    }

    fmpz_init(den_bound);
    fmpz_init(den);

    fmpz_one(den_bound);
    for (i = 0; i <= deg; i++)
    {
        qqbar_denominator(den, coeffs + i);
        fmpz_lcm(den_bound, den_bound, den);
    }
    fmpz_pow_ui(den_bound, den_bound, conjugations);

    /* We have rational coefficients; run the specialized algorithm. */
    if (conjugations == 1)
    {
        fmpz_poly_t t;
        fmpz_poly_init(t);

        fmpz_poly_fit_length(t, deg + 1);
        _fmpz_poly_set_length(t, deg + 1);

        for (i = 0; i <= deg; i++)
        {
            fmpz_divexact(t->coeffs + i, den_bound, QQBAR_POLY(coeffs + i)->coeffs + 1);
            fmpz_mul(t->coeffs + i, t->coeffs + i, QQBAR_POLY(coeffs + i)->coeffs);
            fmpz_neg(t->coeffs + i, t->coeffs + i);
        }

        qqbar_roots_fmpz_poly(roots, t, 0);

        fmpz_poly_clear(t);
        fmpz_clear(den);
        fmpz_clear(den_bound);
        return 1;
    }

    deg_bound = conjugations * deg;
    if (deg_bound > deg_limit)
    {
        fmpz_clear(den);
        fmpz_clear(den_bound);
        return 0;
    }

    conj_roots = flint_malloc(sizeof(acb_ptr) * (deg + 1));

    for (i = 0; i <= deg; i++)
        conj_roots[i] = _acb_vec_init(qqbar_degree(coeffs + i));

    acb_poly_init(cpoly);
    acb_poly_init(ct);
    fmpz_poly_init(rpoly);

    fmpz_one(den_bound);
    for (i = 0; i <= deg; i++)
    {
        qqbar_denominator(den, coeffs + i);
        fmpz_lcm(den_bound, den_bound, den);
    }
    fmpz_pow_ui(den_bound, den_bound, conjugations);

#if VERBOSE
    flint_printf("deg_bound: %wd\n", deg_bound);
    flint_printf("den_bound: %wd bits\n", fmpz_bits(den_bound));
#endif

    initial_prec = 64 + fmpz_bits(den_bound) + 2 * FLINT_BIT_COUNT(conjugations);
    initial_prec = FLINT_MAX(initial_prec, QQBAR_DEFAULT_PREC);

    for (prec = initial_prec; ; prec *= 2)
    {
#if VERBOSE
        flint_printf("prec = %wd\n", prec);
#endif

        for (i = 0; i <= deg; i++)
            arb_fmpz_poly_complex_roots(conj_roots[i], QQBAR_POLY(coeffs + i), 0, prec);

        bsprod(cpoly, 0, conjugations, deg, coeffs, conj_roots, prec);
        _acb_vec_scalar_mul_fmpz(cpoly->coeffs, cpoly->coeffs, cpoly->length, den_bound, prec);

#if VERBOSE
        acb_poly_printd(cpoly, 20); flint_printf("\n\n");
#endif

        if (acb_poly_get_unique_fmpz_poly(rpoly, cpoly))
            break;

        /* todo: reuse enclosures (refine numerical roots) */
        /* todo: compute next precision from coefficient magnitude */


        if (prec > bits_limit)
        {
            success = 0;
            break;
        }
    }

    if (success)
    {
    #if VERBOSE
        flint_printf("rpoly: "); fmpz_poly_print(rpoly); flint_printf("\n\n");
    #endif

        slong found;
        slong num_candidates = fmpz_poly_degree(rpoly);
        slong num_candidates2;

        fmpz_poly_factor_t fac;
        fmpz_poly_factor_init(fac);
        fmpz_poly_factor(fac, rpoly);

        acb_t y;
        acb_init(y);

        qqbar_ptr candidates;
        slong * c_indices;
        slong ci;

        candidates = _qqbar_vec_init(num_candidates);
        c_indices = flint_malloc(sizeof(slong) * num_candidates);

        ci = 0;
        for (i = 0; i < fac->num; i++)
        {
    #if VERBOSE
            flint_printf("factor: ");
            fmpz_poly_print(fac->p + i); flint_printf("\n");
    #endif

            qqbar_roots_fmpz_poly(candidates + ci, fac->p + i, QQBAR_ROOTS_IRREDUCIBLE | QQBAR_ROOTS_UNSORTED);
            ci += fmpz_poly_degree(fac->p + i);
        }

        num_candidates2 = ci;

    #if VERBOSE
        for (ci = 0; ci < num_candidates2; ci++)
        {
            flint_printf("candidate %wd: ", ci);
            qqbar_print(candidates + ci); flint_printf("\n");
        }
    #endif

        initial_prec = QQBAR_DEFAULT_PREC;
        for (prec = initial_prec; ; prec *= 2)
        {
            /* Check numerically for solutions to f(x) = 0. We will necessarily have
               found >= deg. If found == deg, we have certainly isolated the
               correct roots. Otherwise, we retry with higher precision.
               This way we do not need exact qqbar arithmetic to validate f(x) = 0. */
            found = 0;

            acb_poly_fit_length(cpoly, deg + 1);
            _acb_poly_set_length(cpoly, deg + 1);
            for (i = 0; i <= deg; i++)
            {
                /* Todo: copy the original enclosures and reuse for refinement. */
                qqbar_get_acb(cpoly->coeffs + i, coeffs + i, prec);
            }

            for (ci = 0; ci < num_candidates2; ci++)
            {
                /* Refine the candidate to the new precision */
                if (prec != initial_prec)
                    _qqbar_enclosure_raw(QQBAR_ENCLOSURE(candidates + ci), QQBAR_POLY(candidates + ci), QQBAR_ENCLOSURE(candidates + ci), prec);

                acb_poly_evaluate(y, cpoly, QQBAR_ENCLOSURE(candidates + ci), prec);

                if (acb_contains_zero(y))
                {
                    /* -- restart with higher precision */
                    if (found == deg)
                    {
                        found = deg + 1;
                        break;
                    }

                    c_indices[found] = ci;
                    found++;
                }
            }

    #if VERBOSE
            flint_printf("prec %wd : found %wd / deg %wd\n", prec, found, deg);
    #endif

            FLINT_ASSERT(found >= deg);

            if (found == deg)
            {
                for (i = 0; i < deg; i++)
                    qqbar_set(roots + i, candidates + c_indices[i]);

                break;
            }

            if (prec > bits_limit)
            {
                success = 0;
                break;
            }
        }

        _qqbar_vec_clear(candidates, fmpz_poly_degree(rpoly));
        flint_free(c_indices);

        fmpz_poly_factor_clear(fac);
        acb_clear(y);
    }

    /* This could be optional, but it's probably cheap compared to the actual computations... */
    if (success)
        qsort(roots, deg, sizeof(qqbar_struct), (int (*)(const void *, const void *)) qqbar_cmp_root_order);

    for (i = 0; i <= deg; i++)
        _acb_vec_clear(conj_roots[i], qqbar_degree(coeffs + i));

    acb_poly_clear(cpoly);
    acb_poly_clear(ct);
    fmpz_poly_clear(rpoly);
    fmpz_clear(den);
    fmpz_clear(den_bound);
    flint_free(conj_roots);

    return success;
}

