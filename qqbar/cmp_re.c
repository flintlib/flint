/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "qqbar.h"

int
qqbar_cmp_re(const qqbar_t x, const qqbar_t y)
{
    slong prec;
    acb_t z1, z2;
    int res;

    if (!arb_overlaps(acb_realref(QQBAR_ENCLOSURE(x)), acb_realref(QQBAR_ENCLOSURE(y))))
    {
        return arf_cmp(arb_midref(acb_realref(QQBAR_ENCLOSURE(x))),
                       arb_midref(acb_realref(QQBAR_ENCLOSURE(y))));
    }

    if (qqbar_sgn_re(y) == 0)
        return qqbar_sgn_re(x);

    if (qqbar_sgn_re(x) == 0)
        return -qqbar_sgn_re(y);

    if (qqbar_degree(x) == 1 && qqbar_degree(y) == 1)
    {
        /* Reversed order since the signs are reversed */
        return _fmpq_cmp(QQBAR_COEFFS(y), QQBAR_COEFFS(y) + 1,
                         QQBAR_COEFFS(x), QQBAR_COEFFS(x) + 1);
    }

    /* Likely complex conjugates */
    if (fmpz_poly_equal(QQBAR_POLY(x), QQBAR_POLY(y)))
    {
        qqbar_t t;

        /* Complex conjugate quadratics */
        if (qqbar_degree(x) == 2 &&
            !arb_overlaps(acb_imagref(QQBAR_ENCLOSURE(x)), acb_imagref(QQBAR_ENCLOSURE(y))))
            return 0;

        qqbar_init(t);
        qqbar_conj(t, y);
        res = qqbar_equal(x, t);
        qqbar_clear(t);
        if (res == 1)
            return 0;
    }

    /* Subtraction is a scalar operation and will be quick */
    if ((qqbar_degree(x) == 1 || qqbar_degree(y) == 1))
    {
        qqbar_t t;
        qqbar_init(t);
        qqbar_sub(t, x, y);

        res = qqbar_sgn_re(t);

        qqbar_clear(t);
        return res;
    }

    acb_init(z1);
    acb_init(z2);

    acb_set(z1, QQBAR_ENCLOSURE(x));
    acb_set(z2, QQBAR_ENCLOSURE(y));

    res = 0;
    for (prec = QQBAR_DEFAULT_PREC; ; prec *= 2)
    {
        _qqbar_enclosure_raw(z1, QQBAR_POLY(x), z1, prec);
        _qqbar_enclosure_raw(z2, QQBAR_POLY(y), z2, prec);

        if (!arb_overlaps(acb_realref(z1), acb_realref(z2)))
        {
            res = arf_cmp(arb_midref(acb_realref(z1)), arb_midref(acb_realref(z2)));
            break;
        }

        /* Force an exact computation (may be slow) */
        if (prec >= 4 * QQBAR_DEFAULT_PREC)
        {
            qqbar_t t;
            qqbar_init(t);
            qqbar_sub(t, x, y);
            res = qqbar_sgn_re(t);
            qqbar_clear(t);
            break;
        }
    }

    acb_clear(z1);
    acb_clear(z2);

    return res;
}

