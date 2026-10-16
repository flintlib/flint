/*
    Copyright (C) 2026 Edgar Costa

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "fmpq.h"

int
arb_get_simplest_fmpq(fmpq_t res, const arb_t x)
{
    arf_t lo_arf, hi_arf;
    fmpq_t lo, hi;
    int negate;

    if (!arb_is_finite(x))
        return 0;

    /* Fast path: exact ball. The midpoint is the only element of the
     * interval, so it is trivially the simplest rational. arf is
     * dyadic, so arf_get_fmpq produces the canonical form exactly. */
    if (arb_is_exact(x))
    {
        arf_get_fmpq(res, arb_midref(x));
        return 1;
    }

    /* Fast path: 0 lies in the ball. 0/1 has the smallest possible
     * denominator and smallest absolute numerator, so it wins both
     * sub-criteria. arb_contains_zero is an exact predicate on the
     * arf/mag pair, avoiding any fmpq conversion. */
    if (arb_contains_zero(x))
    {
        fmpq_zero(res);
        return 1;
    }

    /* The remaining intervals are strictly positive or strictly
     * negative. Reduce to the positive case so that
     * fmpq_simplest_between's "smallest value" tie-break agrees with
     * our "smallest absolute numerator" contract on the original. */
    negate = arb_is_negative(x);

    arf_init(lo_arf);
    arf_init(hi_arf);
    fmpq_init(lo);
    fmpq_init(hi);

    arb_get_lbound_arf(lo_arf, x, ARF_PREC_EXACT);
    arb_get_ubound_arf(hi_arf, x, ARF_PREC_EXACT);
    arf_get_fmpq(lo, lo_arf);
    arf_get_fmpq(hi, hi_arf);

    if (negate)
    {
        fmpq_neg(lo, lo);
        fmpq_neg(hi, hi);
        fmpq_swap(lo, hi);
    }

    fmpq_simplest_between(res, lo, hi);
    if (negate)
        fmpq_neg(res, res);

    arf_clear(lo_arf);
    arf_clear(hi_arf);
    fmpq_clear(lo);
    fmpq_clear(hi);
    return 1;
}
