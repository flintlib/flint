/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"

slong
acb_rel_error_bits(const acb_t x)
{
    int am, ar, bm, br;
    slong result;
    const fmpz * radmag;
    const fmpz * midmag;
    fmpz_t t;

    am = !arf_is_zero(arb_midref(acb_realref(x)));
    ar = !mag_is_zero(arb_radref(acb_realref(x)));

    bm = !arf_is_zero(arb_midref(acb_imagref(x)));
    br = !mag_is_zero(arb_radref(acb_imagref(x)));

    /* no radius -- exact */
    if (!ar && !br)
        return -ARF_PREC_EXACT;

    /* no midpoint -- infinite relative error */
    if (!am && !bm)
        return ARF_PREC_EXACT;

    if (!acb_is_finite(x))
        return ARF_PREC_EXACT;

#define ame ARF_EXPREF(arb_midref(acb_realref(x)))
#define are MAG_EXPREF(arb_radref(acb_realref(x)))
#define bme ARF_EXPREF(arb_midref(acb_imagref(x)))
#define bre MAG_EXPREF(arb_radref(acb_imagref(x)))

    if (am && bm)
        midmag = fmpz_cmp(ame, bme) >= 0 ? ame : bme;
    else if (am)
        midmag = ame;
    else
        midmag = bme;

    if (ar && br)
        radmag = fmpz_cmp(are, bre) >= 0 ? are : bre;
    else if (ar)
        radmag = are;
    else
        radmag = bre;

    fmpz_init(t);
    fmpz_add_ui(t, radmag, 1);
    result = _fmpz_sub_small(t, midmag);
    fmpz_clear(t);

    return result;
}

slong acb_rel_one_accuracy_bits(const acb_t z)
{
    if (arf_cmpabs_2exp_si(arb_midref(acb_realref(z)), -1) < 0 &&
        arf_cmpabs_2exp_si(arb_midref(acb_imagref(z)), -1) < 0)
    {
        acb_t t;
        arf_init(arb_midref(acb_realref(t)));
        arf_one(arb_midref(acb_realref(t)));
        arf_init(arb_midref(acb_imagref(t)));
        *arb_radref(acb_realref(t)) = *arb_radref(acb_realref(z));
        *arb_radref(acb_imagref(t)) = *arb_radref(acb_imagref(z));
        return acb_rel_accuracy_bits(t);
    }
    else
    {
        return acb_rel_accuracy_bits(z);
    }
}

