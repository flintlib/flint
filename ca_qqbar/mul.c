/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca_qqbar.h"

void
ca_qqbar_mul(ca_qqbar_t res, const ca_qqbar_t x, const ca_qqbar_t y)
{
    if (ca_qqbar_is_zero(x) || ca_qqbar_is_zero(y))
    {
        ca_qqbar_zero(res);
    }
    else if (ca_qqbar_is_one(x))
    {
        ca_qqbar_set(res, y);
    }
    else if (ca_qqbar_is_one(y))
    {
        ca_qqbar_set(res, x);
    }
    else if (ca_qqbar_is_neg_one(x))
    {
        ca_qqbar_neg(res, y);
    }
    else if (ca_qqbar_is_neg_one(y))
    {
        ca_qqbar_neg(res, x);
    }
    else if (ca_qqbar_is_rational(y))
    {
        fmpz_t a, b, c;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);

        fmpz_neg(a, CA_QQBAR_COEFFS(y));
        fmpz_set(c, CA_QQBAR_COEFFS(y) + 1);

        ca_qqbar_scalar_op(res, x, a, b, c);

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
    }
    else if (ca_qqbar_is_rational(x))
    {
        fmpz_t a, b, c;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);

        fmpz_neg(a, CA_QQBAR_COEFFS(x));
        fmpz_set(c, CA_QQBAR_COEFFS(x) + 1);

        ca_qqbar_scalar_op(res, y, a, b, c);

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
    }
    else if (ca_qqbar_equal(x, y))
    {
        /* This may detect exact square roots and other special cases. */
        ca_qqbar_pow_ui(res, x, 2);
    }
    else
    {
        ca_qqbar_binary_op(res, x, y, 2);
    }
}

