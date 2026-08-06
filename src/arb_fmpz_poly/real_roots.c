/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "profiler.h"
#include "fmpz_poly.h"
#include "fmpz_vec.h"
#include "fmpq_vec.h"
#include "acb_poly.h"
#include "arb_fmpz_poly.h"

#ifndef __compar_fn_t
#if defined(_MSC_VER)
typedef int(*__compar_fn_t) (const void *, const void *);
#else
typedef int(*__compar_fn_t) (__const void *, __const void *);
#endif
#endif

static int arb_cmp_mid(const arb_t a, const arb_t b)
{
    return arf_cmp(arb_midref(a), arb_midref(b));
}

static void _arb_vec_sort_mid(arb_ptr vec, slong len)
{
    qsort(vec, len, sizeof(arb_struct), (__compar_fn_t) arb_cmp_mid);
}

static slong
arb_fmpz_poly_real_roots_deflated(arb_ptr res, const fmpz_poly_t poly, slong deflation, int have_zero, int flags, slong prec)
{
    slong deg = fmpz_poly_degree(poly);
    slong * iv_exp, n_ex, n_iv;
    slong i, num_deflated, num_roots;
    fmpz_poly_t q;
    fmpz_t t;
    int positive_only;
    arb_ptr roots_deflated;

    fmpq * ex = _fmpq_vec_init(deg);
    fmpz * iv = _fmpz_vec_init(deg);
    iv_exp = flint_malloc(deg * sizeof(slong));
    fmpz_poly_init(q);
    fmpz_init(t);

    positive_only = (deflation > 1) && (deflation % 2 == 0);

    if (positive_only)
        fmpz_poly_isolate_positive_roots(ex, &n_ex, iv, iv_exp, &n_iv, poly);
    else
        fmpz_poly_isolate_real_roots(ex, &n_ex, iv, iv_exp, &n_iv, poly);

    if (flags & ARB_FMPZ_POLY_ROOTS_VERBOSE)
        flint_printf("isolation gave %wd exact, %wd intervals\n", n_ex, n_iv);

    num_deflated = n_ex + n_iv;
    num_roots = 0;

    roots_deflated = _arb_vec_init(num_deflated);

    for (prec = FLINT_MAX(prec, 30) + ((deflation == 1) ? 0 : 4); ; prec = prec * 2 + 10)
    {
        if (flags & ARB_FMPZ_POLY_ROOTS_VERBOSE)
            flint_printf("Trying precision %wd\n", prec);

        /* Add the exact roots (before inflation) */
        for (i = 0; i < n_ex; i++)
        {
            arb_ptr r = roots_deflated + i;

            //if (flags & ARB_FMPZ_POLY_ROOTS_VERBOSE)
            //    flint_printf("exact %wd: %{fmpq}\n", i, ex + i);

            arb_set_fmpz(r, fmpq_numref(ex + i));
            arb_mul_2exp_si(r, r, -(fmpz_bits(fmpq_denref(ex + i)) - 1));
        }

        /* The output of fmpz_poly_isolate_real_roots is an open interval
           containing exactly one root. However, arb_fmpz_poly_refine_root
           takes a closed interval as input. To inflate (a, b) -> [a, b]
           and refine, we must ensure that a and b are not roots.
           This is most easily done by dividing out possible
           exact roots beforehand.

           XXX: fmpz_poly_isolate_real_roots should typically detect
           exact roots occurring exactly at bisection points, but this is
           not guaranteed by its contract. We may want some extra checks in this
           method to prevent potential hangs.
        */

        if (n_ex != 0)
        {
            fmpz_poly_product_roots_fmpq_vec(q, ex, n_ex);
            fmpz_poly_divexact(q, poly, q);
        }

        /* Refine and add the inexact roots (before inflation) */
        for (i = 0; i < n_iv; i++)
        {
            arb_ptr r = roots_deflated + n_ex + i;

            //if (flags & ARB_FMPZ_POLY_ROOTS_VERBOSE)
            //    flint_printf("interval %wd: %{fmpz} * 2^%wd       prec %wd\n",
            //        i, iv + i, iv_exp[i], prec);

            /* [a, a+1] * 2^k -> [2a+1 +/- 1] * 2^(k-1)  */
            fmpz_mul_2exp(t, iv + i, 1);
            fmpz_add_ui(t, t, 1);
            arb_set_fmpz(r, t);
            mag_one(arb_radref(r));
            arb_mul_2exp_si(r, r, iv_exp[i] - 1);

            arb_fmpz_poly_refine_root_arb(r, (n_ex == 0) ? poly : q, r, prec);
        }

        /* Inflate if needed */
        if (deflation == 1)
        {
            _arb_vec_set(res, roots_deflated, num_deflated);
            num_roots = num_deflated;
        }
        else
        {
            //if (flags & ARB_FMPZ_POLY_ROOTS_VERBOSE)
            //{
            //    flint_printf("deflated roots:\n");
            //    for (i = 0; i < num_deflated; i++)
            //    {
            //        arb_printd(roots_deflated + i, 15);
            //        flint_printf("\n");
            //    }
            //}

            num_roots = 0;

            for (i = 0; i < num_deflated; i++)
            {
                arb_ptr x = roots_deflated + i;

                if (arb_contains_zero(x))
                    flint_throw(FLINT_ERROR, "arb_fmpz_poly_real_roots: root not isolated from 0\n");

                if (!arb_is_exact(x))
                    prec = FLINT_MAX(prec, arb_rel_accuracy_bits(x));

                if (arf_sgn(arb_midref(x)) > 0)
                {
                    arb_root_ui(res + num_roots, x, deflation, prec);
                    num_roots++;
                    if (deflation % 2 == 0)
                    {
                        arb_neg(res + num_roots, res + num_roots - 1);
                        num_roots++;
                    }
                }
                else if (deflation % 2 == 1)
                {
                    arb_neg(res + num_roots, x);
                    arb_root_ui(res + num_roots, res + num_roots, deflation, prec);
                    arb_neg(res + num_roots, res + num_roots);
                    num_roots++;
                }
            }
        }

        if (have_zero)
        {
            arb_zero(res + num_roots);
            num_roots++;
        }

        int isolated = 1;

        _arb_vec_sort_mid(res, num_roots);

        /* It is possible that that the post-deflation operations caused two
           initially isolated intervals to overlap again. We could recover from
           this, but easier will be to just restart from scratch with higher precision. */
        for (i = 0; i < num_roots - 1; i++)
        {
            if (arb_overlaps(res + i, res + i + 1))
            {
                if (flags & ARB_FMPZ_POLY_ROOTS_VERBOSE)
                    flint_printf("UNLIKELY: roots not isolated at precision %wd; restarting\n", prec);

                isolated = 0;
                break;
            }
        }

        if (isolated)
            break;        
    }

    _arb_vec_clear(roots_deflated, num_deflated);

    _fmpq_vec_clear(ex, deg);
    _fmpz_vec_clear(iv, deg);
    flint_free(iv_exp);
    fmpz_poly_clear(q);
    fmpz_clear(t);

    return num_roots;
}


slong
arb_fmpz_poly_real_roots(arb_ptr roots, const fmpz_poly_t poly, int flags, slong target_prec)
{
    slong deflation, num;
    fmpz_poly_t poly_deflated;
    int have_zero;

    if (fmpz_poly_degree(poly) < 1)
        return 0;

    fmpz_poly_init(poly_deflated);

    have_zero = fmpz_is_zero(poly->coeffs);
    if (have_zero)
        fmpz_poly_shift_right(poly_deflated, poly, 1);
    else
        fmpz_poly_set(poly_deflated, poly);

    deflation = arb_fmpz_poly_deflation(poly_deflated);

    arb_fmpz_poly_deflate(poly_deflated, poly_deflated, deflation);

    if (flags & ARB_FMPZ_POLY_ROOTS_VERBOSE)
        flint_printf("degree %wd, %wd deflated\n", fmpz_poly_degree(poly), fmpz_poly_degree(poly_deflated));

    num = arb_fmpz_poly_real_roots_deflated(roots, poly_deflated, deflation, have_zero, flags, target_prec);

    fmpz_poly_clear(poly_deflated);

    return num;
}

