/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "fmpz_poly.h"
#include "acb.h"
#include "qqbar.h"

int
qqbar_set_fmpz_poly_root(qqbar_t res, const fmpz_poly_t poly, const acb_t z, slong max_prec)
{
    acb_t w;
    int ok;

    if (fmpz_poly_degree(poly) < 1 || fmpz_sgn(poly->coeffs + poly->length - 1) < 0)
        return 0;

    if (fmpz_poly_degree(poly) == 1)
    {
        /* (the enclosure may be arbitrarily wide for a rational root) */
        fmpq_t q;
        fmpq_init(q);
        fmpz_neg(fmpq_numref(q), poly->coeffs);
        fmpz_set(fmpq_denref(q), poly->coeffs + 1);
        fmpq_canonicalise(q);
        ok = acb_contains_fmpq(z, q);
        if (ok)
            qqbar_set_fmpq(res, q);
        fmpq_clear(q);
        return ok;
    }

    acb_init(w);
    ok = _qqbar_validate_existence_uniqueness(w, poly, z, max_prec);

    if (ok)
    {
        fmpz_poly_set(QQBAR_POLY(res), poly);
        acb_set(QQBAR_ENCLOSURE(res), w);
    }

    acb_clear(w);
    return ok;
}
