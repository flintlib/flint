/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "nmod.h"

ulong
nmod_pow_fmpz(ulong a, const fmpz_t exp, nmod_t mod)
{
    FLINT_ASSERT(a < mod.n);

    if (fmpz_sgn(exp) < 0)
    {
        ulong g, r;
        fmpz_t e;

        /* as n_powmod2_preinv: a power of the inverse */
        g = n_gcdinv(&a, a, mod.n);

        if (g != 1)
            flint_throw(FLINT_IMPINV, "Cannot invert modulo %wu*%wu\n",
                        g, mod.n / g);

        if (!COEFF_IS_MPZ(*exp))
            return nmod_pow_ui(a, -(ulong) *exp, mod);

        fmpz_init(e);
        fmpz_neg(e, exp);
        r = n_powmod2_fmpz_preinv(a, e, mod.n, mod.ninv);
        fmpz_clear(e);
        return r;
    }

    return n_powmod2_fmpz_preinv(a, exp, mod.n, mod.ninv);
}
