/*
    Copyright (C) 2015 Kushagra Singh

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "mpn_extras.h"

void
fmpz_factor_ecm_init(ecm_t ecm_inf, mp_limb_t sz)
{
    ecm_inf->t = flint_malloc(sz * sizeof(mp_limb_t));
    ecm_inf->u = flint_malloc(sz * sizeof(mp_limb_t));
    ecm_inf->v = flint_malloc(sz * sizeof(mp_limb_t));
    ecm_inf->w = flint_malloc(sz * sizeof(mp_limb_t));

    ecm_inf->x = flint_malloc(sz * sizeof(mp_limb_t));
    ecm_inf->z = flint_malloc(sz * sizeof(mp_limb_t));

    ecm_inf->a24 = flint_malloc(sz * sizeof(mp_limb_t));
    ecm_inf->ninv = flint_malloc(sz * sizeof(mp_limb_t));
    ecm_inf->one = flint_malloc(sz * sizeof(mp_limb_t));

    mpn_zero(ecm_inf->t, sz);
    mpn_zero(ecm_inf->u, sz);
    mpn_zero(ecm_inf->v, sz);
    mpn_zero(ecm_inf->w, sz);

    mpn_zero(ecm_inf->x, sz);
    mpn_zero(ecm_inf->z, sz);

    mpn_zero(ecm_inf->a24, sz);
    mpn_zero(ecm_inf->ninv, sz);
    mpn_zero(ecm_inf->one, sz);

    ecm_inf->n_size = sz;
}
