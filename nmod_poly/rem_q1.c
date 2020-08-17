/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "ulong_extras.h"

void _nmod_poly_rem_q1(mp_ptr R, 
                       mp_srcptr A, slong lenA, mp_srcptr B, slong lenB,
                       nmod_t mod)
{
    const mp_limb_t invL = (B[lenB-1] == 1) ? 1 : n_invmod(B[lenB-1], mod.n);

    if (lenB > 1)
    {
        mp_limb_t t, q0, q1;

        q1 = n_mulmod2_preinv(A[lenA-1], invL, mod.n, mod.ninv);
        t  = n_mulmod2_preinv(q1, B[lenB-2], mod.n, mod.ninv);
        t  = n_submod(A[lenA-2], t, mod.n);
        q0 = n_mulmod2_preinv(t, invL, mod.n, mod.ninv);

        if (FLINT_BITS + 2 <= 2 * mod.norm)
        {
            mpn_mul_1(R, B, lenB - 1, q0);
            if (lenB > 2)
                mpn_addmul_1(R + 1, B, lenB - 2, q1);
            _nmod_vec_reduce(R, R, lenB - 1, mod);
        }
        else
        {
            _nmod_vec_scalar_mul_nmod(R, B, lenB - 1, q0, mod);
            if (lenB > 2)
                _nmod_vec_scalar_addmul_nmod(R + 1, B, lenB - 2, q1, mod);
        }

        _nmod_vec_sub(R, A, R, lenB - 1, mod);
    }
}

