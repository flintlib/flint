/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arith.h"

void
arith_bell_number_nmod_vec_series(mp_ptr res, slong n, nmod_t mod)
{
    mp_limb_t fac, c;
    mp_ptr tmp;
    slong k;

    if (n < 1)
        return;

    tmp = flint_malloc(sizeof(mp_limb_t) * n);

    /* Divide by factorials */
    fac = n_factorial_mod2_preinv(n-1, mod.n, mod.ninv);
    c = n_invmod(fac, mod.n);

    for (k = n - 1; k > 0; k--)
    {
        tmp[k] = c;
        c = n_mulmod2_preinv(c, k, mod.n, mod.ninv);
    }
    tmp[0] = UWORD(0);

    _nmod_poly_exp_series(res, tmp, n, n, mod);

    /* Multiply by factorials */
    c = UWORD(1);
    for (k = 1; k < n; k++)
    {
        c = n_mulmod2_preinv(c, k, mod.n, mod.ninv);
        res[k] = n_mulmod2_preinv(res[k], c, mod.n, mod.ninv);
    }

    flint_free(tmp);
}
