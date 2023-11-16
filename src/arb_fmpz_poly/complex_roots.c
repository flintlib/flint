/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"
#include "arb_fmpz_poly.h"
#include "profiler.h"

static int check_accuracy(acb_ptr vec, slong len, slong prec)
{
    slong i;

    for (i = 0; i < len; i++)
    {
        if (acb_rel_accuracy_bits(vec + i) < prec)
            return 0;
    }

    return 1;
}

static int check_isolation(acb_srcptr roots, slong len)
{
    slong i, j;

    for (i = 0; i < len; i++)
    {
        if (arf_sgn(arb_midref(acb_imagref(roots + i))) >= 0)
        {
            for (j = i + 1; j < len; j++)
            {
                if (arf_sgn(arb_midref(acb_imagref(roots + j))) >= 0)
                {
                    if (acb_overlaps(roots + i, roots + j))
                        return 0;
                }
            }
        }
    }

    return 1;
}

void
arb_fmpz_poly_complex_roots(acb_ptr roots, const fmpz_poly_t poly, int flags, slong target_prec)
{
    slong i, j, prec, deg, deg_deflated, isolated, maxiter, deflation;
    slong initial_prec, num_real;
    acb_poly_t cpoly, cpoly_deflated;
    fmpz_poly_t poly_deflated;
    acb_ptr roots_deflated;
    int removed_zero;

    if (fmpz_poly_degree(poly) < 1)
        return;

    initial_prec = 32;

    fmpz_poly_init(poly_deflated);
    acb_poly_init(cpoly);
    acb_poly_init(cpoly_deflated);

    /* try to write poly as poly_deflated(x^deflation), possibly multiplied by x */
    removed_zero = fmpz_is_zero(poly->coeffs);
    if (removed_zero)
        fmpz_poly_shift_right(poly_deflated, poly, 1);
    else
        fmpz_poly_set(poly_deflated, poly);
    deflation = arb_fmpz_poly_deflation(poly_deflated);
    arb_fmpz_poly_deflate(poly_deflated, poly_deflated, deflation);

    deg = fmpz_poly_degree(poly);
    deg_deflated = fmpz_poly_degree(poly_deflated);

    if (flags & ARB_FMPZ_POLY_ROOTS_VERBOSE)
    {
        flint_printf("searching for %wd roots, %wd deflated\n", deg, deg_deflated);
    }

    /* we only need deg_deflated entries, but the remainder will be useful
       as scratch space */
    roots_deflated = _acb_vec_init(deg);

    for (prec = initial_prec; ; prec *= 2)
    {
        acb_poly_set_fmpz_poly(cpoly_deflated, poly_deflated, prec);
        maxiter = FLINT_MIN(4 * deg_deflated + 64, prec);

        if (flags & ARB_FMPZ_POLY_ROOTS_VERBOSE)
        {
            TIMEIT_ONCE_START
            flint_printf("prec=%wd: ", prec);
            isolated = acb_poly_find_roots(roots_deflated, cpoly_deflated,
                prec == initial_prec ? NULL : roots_deflated, maxiter, prec);
            flint_printf("%wd isolated roots | ", isolated);
            TIMEIT_ONCE_STOP
        }
        else
        {
            isolated = acb_poly_find_roots(roots_deflated, cpoly_deflated,
                prec == initial_prec ? NULL : roots_deflated, maxiter, prec);
        }

        if (isolated == deg_deflated)
        {
            if (!check_accuracy(roots_deflated, deg_deflated, target_prec))
                continue;

            if (deflation == 1)
            {
                _acb_vec_set(roots, roots_deflated, deg_deflated);
            }
            else  /* compute all nth roots */
            {
                acb_t w, w2;

                acb_init(w);
                acb_init(w2);

                acb_unit_root(w, deflation, prec);
                acb_unit_root(w2, 2 * deflation, prec);

                for (i = 0; i < deg_deflated; i++)
                {
                    if (arf_sgn(arb_midref(acb_realref(roots_deflated + i))) > 0)
                    {
                        acb_root_ui(roots + i * deflation,
                                    roots_deflated + i, deflation, prec);
                    }
                    else
                    {
                        acb_neg(roots + i * deflation, roots_deflated + i);
                        acb_root_ui(roots + i * deflation,
                            roots + i * deflation, deflation, prec);
                        acb_mul(roots + i * deflation,
                            roots + i * deflation, w2, prec);
                    }

                    for (j = 1; j < deflation; j++)
                    {
                        acb_mul(roots + i * deflation + j,
                                roots + i * deflation + j - 1, w, prec);
                    }
                }

                acb_clear(w);
                acb_clear(w2);
            }

            /* by assumption that poly is squarefree, must be just one */
            if (removed_zero)
                acb_zero(roots + deg_deflated * deflation);

            if (!check_accuracy(roots, deg, target_prec))
                continue;

            acb_poly_set_fmpz_poly(cpoly, poly, prec);

            if (!acb_poly_validate_real_roots(roots, cpoly, prec))
                continue;

            for (i = 0; i < deg; i++)
            {
                if (arb_contains_zero(acb_imagref(roots + i)))
                    arb_zero(acb_imagref(roots + i));
            }

            if (!check_isolation(roots, deg))
            {
                /* extremely unlikely */
                if (flags & ARB_FMPZ_POLY_ROOTS_VERBOSE)
                    flint_printf("isolation failure!\n");

                continue;
            }

            if (flags & ARB_FMPZ_POLY_ROOTS_VERBOSE)
                flint_printf("done!\n");

            break;
        }
    }

    _acb_vec_sort_pretty(roots, deg);

    /* pair conjugates */
    num_real = 0;
    for (i = 0; i < deg; i++)
        if (acb_is_real(roots + i))
            num_real++;

    if (deg != num_real)
    {
        for (i = num_real, j = 0; i < deg; i++)
        {
            if (arb_is_positive(acb_imagref(roots + i)))
            {
                acb_swap(roots_deflated + j, roots + i);
                j++;
            }
        }

        for (i = 0; i < (deg - num_real) / 2; i++)
        {
            acb_swap(roots + num_real + 2 * i, roots_deflated + i);
            acb_conj(roots + num_real + 2 * i + 1, roots + num_real + 2 * i);
        }
    }

    fmpz_poly_clear(poly_deflated);
    acb_poly_clear(cpoly);
    acb_poly_clear(cpoly_deflated);
    _acb_vec_clear(roots_deflated, deg);
}

