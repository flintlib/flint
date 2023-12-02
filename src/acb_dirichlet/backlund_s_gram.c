/*
    Copyright (C) 2019 D.H.J. Polymath

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

slong
acb_dirichlet_backlund_s_gram(const fmpz_t n)
{
    slong res = 0;
    if (fmpz_cmp_si(n, -1) < 0)
    {
        flint_throw(FLINT_ERROR, "n must be >= -1\n");
    }
    else
    {
        fmpz_t k;
        fmpz_init(k);
        acb_dirichlet_zeta_nzeros_gram(k, n);
        fmpz_sub(k, k, n);
        res = fmpz_get_si(k) - 1;
        fmpz_clear(k);
    }
    return res;
}
