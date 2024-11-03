/*
    Copyright (C) 2011, 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "arith.h"
#include "arith-impl.h"

ulong nmod_inv_check(ulong x, nmod_t mod)
{
    ulong r, g;

    g = n_gcdinv(&r, x, mod.n);
    if (g != 1)
        return mod.n;

    return r;
}

int
arith_bell_number_nmod_vec_series(nn_ptr res, slong n, nmod_t mod)
{
    ulong c;
    nn_ptr tmp;
    slong k;
    int success;

    if (n <= 0)
        return 1;

    if (mod.n == 1)
        return 0;

    tmp = flint_malloc(sizeof(ulong) * n);

    /* Compute inverse factorials */
    c = 1;
    if (NMOD_BITS(mod) == FLINT_BITS)
    {
        for (k = n - 1; k > 0; k--)
        {
            tmp[k] = c;
            c = _nmod_mul_fullword(c, k, mod);
        }
    }
    else
    {
        for (k = n - 1; k > 0; k--)
        {
            tmp[k] = c;
            c = nmod_mul(c, k, mod);
        }
    }

    c = nmod_inv_check(c, mod);
    success = (c != mod.n);

    if (success)
    {
        tmp[0] = 0;
        _nmod_vec_scalar_mul_nmod(tmp + 1, tmp + 1, n - 1, c, mod);

        _nmod_poly_exp_series(res, tmp, n, n, mod);

        /* Multiply by factorials */
        c = 1;
        if (NMOD_BITS(mod) == FLINT_BITS)
        {
            for (k = 1; k < n; k++)
            {
                c = _nmod_mul_fullword(c, k, mod);
                res[k] = _nmod_mul_fullword(res[k], c, mod);
            }
        }
        else
        {
            for (k = 1; k < n; k++)
            {
                c = nmod_mul(c, k, mod);
                res[k] = nmod_mul(res[k], c, mod);
            }
        }
    }

    flint_free(tmp);
    return success;
}
