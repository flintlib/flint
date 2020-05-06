/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca_qqbar.h"

int
ca_qqbar_acot_pi(slong * p, ulong * q, const ca_qqbar_t x)
{
    slong deg = ca_qqbar_degree(x);

    *p = 0;
    *q = 1;

    if (deg == 1)
    {
        if (ca_qqbar_is_zero(x))
        {
            *p = 1;
            *q = 2;
            return 1;
        }

        if (ca_qqbar_is_one(x))
        {
            *p = 1;
            *q = 4;
            return 1;
        }

        if (ca_qqbar_is_neg_one(x))
        {
            *p = -1;
            *q = 4;
            return 1;
        }

        return 0;
    }
    else if (deg == 2)
    {
        fmpz a, b, c;

        a = CA_QQBAR_COEFFS(x)[0];
        b = CA_QQBAR_COEFFS(x)[1];
        c = CA_QQBAR_COEFFS(x)[2];

        if (a == -3 && b == 0 && c == 1)
        {
            *p = ca_qqbar_sgn_re(x);
            *q = 6;
            return 1;
        }

        if (a == -1 && b == 0 && c == 3)
        {
            *p = ca_qqbar_sgn_re(x);
            *q = 3;
            return 1;
        }

        if (a == -1 && b == 2 && c == 1)
        {
            *p = (ca_qqbar_sgn_re(x) == 1) ? 3 : -1;
            *q = 8;
            return 1;
        }

        if (a == -1 && b == -2 && c == 1)
        {
            *p = (ca_qqbar_sgn_re(x) == 1) ? 1 : -3;
            *q = 8;
            return 1;
        }

        if (a == 1 && b == -4 && c == 1)
        {
            /* root is ~0.267 or ~3.73 -- accuracy should not be that bad */
            if (arb_contains_si(acb_realref(CA_QQBAR_ENCLOSURE(x)), 1))
                flint_abort();
 
           *p = (arf_cmpabs_2exp_si(arb_midref(acb_realref(CA_QQBAR_ENCLOSURE(x))), 0) < 0) ? 5 : 1;
            *q = 12;
            return 1;
        }

        if (a == 1 && b == 4 && c == 1)
        {
            if (arb_contains_si(acb_realref(CA_QQBAR_ENCLOSURE(x)), -1))
                flint_abort();
            *p = (arf_cmpabs_2exp_si(arb_midref(acb_realref(CA_QQBAR_ENCLOSURE(x))), 0) < 0) ? -5 : -1;
            *q = 12;
            return 1;
        }

        return 0;
    }
    else if ((deg % 2 != 0) || !ca_qqbar_is_real(x))
    {
        return 0;
    }
    else
    {
        int res;
        ca_qqbar_t t;
        ca_qqbar_init(t);
        ca_qqbar_inv(t, x);
        res = ca_qqbar_atan_pi(p, q, t);
        ca_qqbar_clear(t);
        return res;
    }
}

