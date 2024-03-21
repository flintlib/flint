/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "ulong_extras.h"

mp_limb_t n_primitive_root_prime_prefactor(ulong p, n_factor_t * factors)
{
    if (p == 2)
        return 1;

    // compute the divisions "(p-1) / factors" once for all
    ulong * exps = FLINT_ARRAY_ALLOC(factors->num, ulong);
    for (slong i = 0; i < factors->num; i++)
        exps[i] = (p-1) / factors->p[i];

    // try 2, 3, ..., p-2
    const ulong pinv = n_preinvert_limb(p);
    for (ulong a = 2; a < p-1; a++)
    {
        slong i = 0;
        while ((i < factors->num) && (1 == n_powmod2_preinv(a, exps[i], p, pinv)))
            i += 1;
        if (i == factors->num)
        {
            flint_free(exps);
            return a;
        }
    }

    flint_free(exps);
    flint_throw(FLINT_ERROR, "Exception (n_primitive_root_prime_prefactor).  root not found.\n");
}

mp_limb_t n_primitive_root_prime(ulong p)
{
    n_factor_t factors;
    n_factor_init(&factors);
    n_factor(&factors, p - 1, 1);

    mp_limb_t a = n_primitive_root_prime_prefactor(p, &factors);
    return a;
}
