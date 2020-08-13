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
fmpz_factor_ecm_clear(ecm_t ecm_inf)
{
    flint_free(ecm_inf->t);
    flint_free(ecm_inf->u);
    flint_free(ecm_inf->v);
    flint_free(ecm_inf->w);

    flint_free(ecm_inf->x);
    flint_free(ecm_inf->z);

    flint_free(ecm_inf->a24);
    flint_free(ecm_inf->ninv);
    flint_free(ecm_inf->one);
}
