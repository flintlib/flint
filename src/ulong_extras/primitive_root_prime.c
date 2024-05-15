/*
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2024 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"

ulong n_primitive_root_prime_prefactor(ulong p, n_factor_t * factors)
{
    if (p == 2)
        return 1;

    // compute the divisions "(p-1) / factors" once for all
    slong exps[FLINT_MAX_FACTORS_IN_LIMB];
    for (slong i = 0; i < factors->num; i++)
        exps[i] = (p-1) / factors->p[i];

    // try 2, 3, ..., p-1
    const ulong pinv = n_preinvert_limb(p);
    for (ulong a = 2; a < p; a++)
    {
        slong i = 0;
        while ((i < factors->num) && (1 != n_powmod2_preinv(a, exps[i], p, pinv)))
            i += 1;
        if (i == factors->num)
            return a;
    }

    flint_throw(FLINT_ERROR, "Exception (n_primitive_root_prime_prefactor).  root not found.\n");
}

ulong n_primitive_root_prime(ulong p)
{
    n_factor_t factors;
    n_factor_init(&factors);
    n_factor(&factors, p - 1, 1);
    return n_primitive_root_prime_prefactor(p, &factors);
}
