/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_fmpz_poly.h"
#include "qqbar.h"

void
_qqbar_enclosure_raw(acb_t res, const fmpz_poly_t poly, const acb_t zin, slong prec)
{
    slong d, orig_prec, step, acc;
    fmpz_poly_t deriv;
    acb_t z, zmid, t, u;

    if (acb_rel_accuracy_bits(zin) >= prec - 3)
    {
        acb_set(res, zin);
        return;
    }

    d = fmpz_poly_degree(poly);

    if (d == 1)
    {
        arb_set_fmpz(acb_realref(res), poly->coeffs);
        arb_div_fmpz(acb_realref(res), acb_realref(res), poly->coeffs + 1, prec);
        arb_neg(acb_realref(res), acb_realref(res));
        arb_zero(acb_imagref(res));
        return;
    }

    orig_prec = prec;

    acc = acb_rel_accuracy_bits(zin);
    prec = FLINT_MAX(acc, 32) + 10;

    fmpz_poly_init(deriv);
    fmpz_poly_derivative(deriv, poly);
    acb_init(z);
    acb_init(zmid);
    acb_init(t);
    acb_init(u);

    acb_set(z, zin);

    for (step = 0; ; step++)
    {
/*
        printf("trying %ld, %ld\n", step, prec);
        acb_printn(z, prec, ARB_STR_CONDENSE * 10);
        printf("\n");
*/

        if (step > 40 || prec > 1000000000)
        {
            flint_throw(FLINT_ERROR, "qqbar_enclosure_raw: root refinement not converging\n");
        }

        prec *= 2;

        /* Interval Newton refinement */
        acb_get_mid(zmid, z);
        arb_fmpz_poly_evaluate_acb(t, poly, zmid, prec);
        arb_fmpz_poly_evaluate_acb(u, deriv, z, prec);
        acb_div(t, t, u, prec);
        acb_sub(t, zmid, t, prec);

        if (acb_contains(z, t) && acb_rel_accuracy_bits(t) >= 1.1 * orig_prec)
        {
            acb_set(res, t);
            break;
        }

        if (acb_contains(z, t) && acb_rel_accuracy_bits(t) > 1.5 * acb_rel_accuracy_bits(z))
        {
            /* Use refined value for next iteration */
            acb_set(z, t);
        }
        else
        {
            /* Newton refinement seems to be converging too slowly --
               force a slow recomputation of all roots (this could be improved...) */
            acb_ptr roots;
            slong found, i;

/*
            printf("recompute %ld %ld\n", d, prec);
*/
            roots = _acb_vec_init(d);

/*
            fmpz_poly_print(poly);
            printf("\n");
*/

            if (!fmpz_poly_is_squarefree(poly))
                flint_throw(FLINT_ERROR, "(%s)\n", __func__);

            arb_fmpz_poly_complex_roots(roots, poly, 0, 2 * prec);

            /* Check for unique root */
            found = -1;
            for (i = 0; i < d && found != -2; i++)
            {
                if (acb_overlaps(roots + i, z))
                {
                    if (found == -1)
                        found = i;
                    else
                        found = -2;
                }
            }

            if (found >= 0)
                acb_set(t, roots + found);

            _acb_vec_clear(roots, d);

            if (found >= 0 && acb_rel_accuracy_bits(t) >= 1.1 * orig_prec)
            {
                acb_set(res, t);
                break;
            }

            if (found >= 0 && acb_rel_accuracy_bits(t) > 1.5 * acb_rel_accuracy_bits(z) + 1)
            {
                /* Use refined value for next iteration */
                acb_set(z, t);
            }
        }
    }

    fmpz_poly_clear(deriv);
    acb_clear(z);
    acb_clear(zmid);
    acb_clear(t);
    acb_clear(u);
}

void
qqbar_enclosure_raw(acb_t res, const qqbar_t x, slong prec)
{
    _qqbar_enclosure_raw(res, QQBAR_POLY(x), QQBAR_ENCLOSURE(x), prec);
}

