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
ca_qqbar_sin_pi(ca_qqbar_t res, slong p, ulong q)
{
    /* we could allow larger unreduced input, but it's a hassle... */
    if (FLINT_BIT_COUNT(q) > FLINT_BITS - 4 || FLINT_BIT_COUNT(FLINT_ABS(p)) > FLINT_BITS - 4)
    {
        flint_printf("ca_qqbar_sin_pi: unexpected large input\n");
        flint_abort();
    }

    p = p % (2 * (slong) q);
    ca_qqbar_cos_pi(res, (slong) q - 2 * p, 2 * q);
}

