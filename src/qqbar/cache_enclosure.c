/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "qqbar.h"

void
qqbar_cache_enclosure(qqbar_t res, slong prec)
{
    acb_t t;
    slong want_prec;

    want_prec = FLINT_MAX(QQBAR_DEFAULT_PREC, prec) * 1.1 + 32;

    acb_init(t);
    qqbar_get_acb(t, res, want_prec);

    if (acb_contains(QQBAR_ENCLOSURE(res), t))
        acb_swap(QQBAR_ENCLOSURE(res), t);

    acb_clear(t);
}

