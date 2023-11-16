/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arf.h"

int
arf_equal(const arf_t x, const arf_t y)
{
    mp_size_t n;

    if (x == y)
        return 1;

    if (ARF_XSIZE(x) != ARF_XSIZE(y))
        return 0;

    if (!fmpz_equal(ARF_EXPREF(x), ARF_EXPREF(y)))
        return 0;

    n = ARF_SIZE(x);

    if (n == 0)
        return 1;

    if (n == 1)
        return ARF_NOPTR_D(x)[0] == ARF_NOPTR_D(y)[0];

    if (n == 2)
        return (ARF_NOPTR_D(x)[0] == ARF_NOPTR_D(y)[0] &&
                ARF_NOPTR_D(x)[1] == ARF_NOPTR_D(y)[1]);

    return mpn_cmp(ARF_PTR_D(x), ARF_PTR_D(y), n) == 0;
}

int
arf_equal_si(const arf_t x, slong y)
{
    arf_t t;
    arf_init_set_si(t, y);
    return arf_equal(x, t); /* no need to free */
}

int
arf_equal_ui(const arf_t x, ulong y)
{
    arf_t t;
    arf_init_set_ui(t, y);
    return arf_equal(x, t); /* no need to free */
}

int
arf_equal_d(const arf_t x, double y)
{
    arf_t t;
    arf_init(t);
    arf_set_d(t, y);
    return arf_equal(x, t); /* no need to free */
}
