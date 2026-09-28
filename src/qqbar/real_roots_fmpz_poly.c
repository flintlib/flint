/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "fmpq.h"
#include "fmpz_poly.h"
#include "fmpz_poly_factor.h"
#include "fmpq_poly.h"
#include "arb_fmpz_poly.h"
#include "qqbar.h"

/* poly irreducible of degree >= 2 */
static slong
_qqbar_real_roots_fmpz_poly_irreducible(qqbar_ptr res, const fmpz_poly_t poly)
{
    slong d, r, i, prec;
    fmpz_t c;
    arb_ptr roots;
    acb_t t;

    d = fmpz_poly_degree(poly);

    /* Odd degree polynomials always have a real root; for even degree,
       we could check quickly for no real roots, but real root isolation
       handles this efficiently. */
    roots = _arb_vec_init(d);
    acb_init(t);
    fmpz_init(c);

    fmpz_poly_content(c, poly);
    if (fmpz_sgn(poly->coeffs + d) < 0)
        fmpz_neg(c, c);

    for (prec = QQBAR_DEFAULT_PREC; ; prec *= 2)
    {
        r = arb_fmpz_poly_real_roots(roots, poly, 0, prec);

        for (i = 0; i < r; i++)
        {
            acb_set_arb(t, roots + i);
            if (!_qqbar_validate_uniqueness(t, poly, t, prec))
                break;
            arb_set(roots + i, acb_realref(t));
        }

        if (i == r)
            break;
    }

    for (i = 0; i < r; i++)
    {
        if (fmpz_is_one(c))
            fmpz_poly_set(QQBAR_POLY(res + i), poly);
        else
            fmpz_poly_scalar_divexact_fmpz(QQBAR_POLY(res + i), poly, c);

        acb_set_arb(QQBAR_ENCLOSURE(res + i), roots + i);
    }

    _arb_vec_clear(roots, d);
    acb_clear(t);
    fmpz_clear(c);

    return r;
}

static void
_qqbar_set_linear_root(qqbar_t res, const fmpz_poly_t poly)
{
    fmpq_t t;
    fmpq_init(t);
    fmpz_neg(fmpq_numref(t), poly->coeffs);
    fmpz_set(fmpq_denref(t), poly->coeffs + 1);
    fmpq_canonicalise(t);
    qqbar_set_fmpq(res, t);
    fmpq_clear(t);
}

slong
qqbar_real_roots_fmpz_poly(qqbar_ptr res, const fmpz_poly_t poly, int flags)
{
    slong d, r;

    d = fmpz_poly_degree(poly);

    if (d <= 0)
        return 0;

    if (d == 1)
    {
        _qqbar_set_linear_root(res, poly);
        return 1;
    }

    if (flags & QQBAR_ROOTS_IRREDUCIBLE)
    {
        r = _qqbar_real_roots_fmpz_poly_irreducible(res, poly);
    }
    else
    {
        fmpz_poly_factor_t sqf, fac;
        slong i, j, k, e, rj;
        qqbar_ptr out;

        fmpz_poly_factor_init(sqf);
        fmpz_poly_factor_init(fac);

        /* Only fully factor the squarefree parts which have real roots. */
        fmpz_poly_factor_squarefree(sqf, poly);

        r = 0;
        out = res;

        for (i = 0; i < sqf->num; i++)
        {
            if (fmpz_poly_degree(sqf->p + i) >= 2 &&
                fmpz_poly_num_real_roots(sqf->p + i) == 0)
                continue;

            e = sqf->exp[i];

            fmpz_poly_factor(fac, sqf->p + i);

            for (j = 0; j < fac->num; j++)
            {
                if (fmpz_poly_degree(fac->p + j) == 1)
                {
                    _qqbar_set_linear_root(out, fac->p + j);
                    rj = 1;
                }
                else
                {
                    rj = _qqbar_real_roots_fmpz_poly_irreducible(out, fac->p + j);
                }

                /* duplicate entries with higher multiplicity */
                if (e > 1)
                {
                    for (k = rj * e - 1; k > 0; k--)
                        qqbar_set(out + k, out + k / e);
                }

                out += rj * e;
                r += rj * e;
            }
        }

        fmpz_poly_factor_clear(sqf);
        fmpz_poly_factor_clear(fac);
    }

    if (!(flags & QQBAR_ROOTS_UNSORTED) && r > 1)
        qsort(res, r, sizeof(qqbar_struct), (int (*)(const void *, const void *)) qqbar_cmp_root_order);

    return r;
}

slong
qqbar_real_roots_fmpq_poly(qqbar_ptr res, const fmpq_poly_t poly, int flags)
{
    fmpz_poly_t t;
    slong r;

    t->coeffs = poly->coeffs;
    t->length = poly->length;
    t->alloc = poly->alloc;

    r = qqbar_real_roots_fmpz_poly(res, t, flags);

    return r;
}
