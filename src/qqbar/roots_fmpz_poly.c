/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "fmpq.h"
#include "fmpz_poly_factor.h"
#include "arb_fmpz_poly.h"
#include "qqbar.h"

void
qqbar_roots_fmpz_poly(qqbar_ptr res, const fmpz_poly_t poly, int flags)
{
    slong d = fmpz_poly_degree(poly);

    if (d == 0 || d == -1)
        return;

    if (d == 1)
    {
        fmpq_t t;
        fmpq_init(t);
        fmpz_neg(fmpq_numref(t), poly->coeffs);
        fmpz_set(fmpq_denref(t), poly->coeffs + 1);
        fmpq_canonicalise(t);  /* irreducible, but may still have content */
        qqbar_set_fmpq(res, t);
        fmpq_clear(t);
        return;
    }

    if (flags & QQBAR_ROOTS_IRREDUCIBLE)
    {
        slong prec, i, checked;
        fmpz_t c;
        acb_ptr croots;

        croots = _acb_vec_init(d);
        fmpz_init(c);
        fmpz_poly_content(c, poly);
        if (fmpz_sgn(poly->coeffs + d) < 0)
            fmpz_neg(c, c);

        for (prec = QQBAR_DEFAULT_PREC; ; prec *= 2)
        {
            arb_fmpz_poly_complex_roots(croots, poly, 0, prec);

            checked = 0;
            for (i = 0; i < d; i++)
            {
                if (_qqbar_validate_uniqueness(croots + i, poly, croots + i, prec))
                    checked++;
                else
                    break;
            }

            if (checked == d)
            {
                for (i = 0; i < d; i++)
                {
                    if (fmpz_is_one(c))
                        fmpz_poly_set(QQBAR_POLY(res + i), poly);
                    else
                        fmpz_poly_scalar_divexact_fmpz(QQBAR_POLY(res + i), poly, c);
                    acb_set(QQBAR_ENCLOSURE(res + i), croots + i);
                }

                break;
            }
        }

        _acb_vec_clear(croots, d);
        fmpz_clear(c);
    }
    else
    {
        fmpz_poly_factor_t fac;
        qqbar_ptr out;
        slong i, j, k, e, facd;

        fmpz_poly_factor_init(fac);
        fmpz_poly_factor(fac, poly);

        out = res;
        for (i = 0; i < fac->num; i++)
        {
            facd = fmpz_poly_degree(fac->p + i);
            qqbar_roots_fmpz_poly(out, fac->p + i, QQBAR_ROOTS_IRREDUCIBLE);
            e = fac->exp[i];

            /* duplicate entries with higher multiplicity */
            if (e > 1)
            {
                for (j = facd - 1; j >= 0; j--)
                {
                    qqbar_set(out + j * e, out + j);
                    for (k = 1; k < e; k++)
                        qqbar_set(out + j * e + k, out + j * e);
                }
            }

            out += e * facd;
        }

        fmpz_poly_factor_clear(fac);
    }

    if (!(flags & QQBAR_ROOTS_UNSORTED))
    {
        qsort(res, d, sizeof(qqbar_struct), (int (*)(const void *, const void *)) qqbar_cmp_root_order);
    }
}
