/*
    Copyright (C) 2015 Kushagra Singh

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_factor.h"

void
fmpz_factor_ecm_init(ecm_t ecm_inf, ulong sz)
{
    ecm_inf->t = flint_calloc(sz, sizeof(ulong));
    ecm_inf->u = flint_calloc(sz, sizeof(ulong));
    ecm_inf->v = flint_calloc(sz, sizeof(ulong));
    ecm_inf->w = flint_calloc(sz, sizeof(ulong));

    ecm_inf->x = flint_calloc(sz, sizeof(ulong));
    ecm_inf->z = flint_calloc(sz, sizeof(ulong));

    ecm_inf->a24 = flint_calloc(sz, sizeof(ulong));
    ecm_inf->ninv = flint_calloc(sz, sizeof(ulong));
    ecm_inf->one = flint_calloc(sz, sizeof(ulong));

    ecm_inf->n_size = sz;
}
