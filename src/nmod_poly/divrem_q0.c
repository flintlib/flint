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

void _nmod_poly_divrem_q0(mp_ptr Q, mp_ptr R, 
                          mp_srcptr A, mp_srcptr B, slong lenA, nmod_t mod)
{
    const mp_limb_t invL = (B[lenA-1] == 1) ? 1 : n_invmod(B[lenA-1], mod.n);

    if (lenA == 1)
    {
        _nmod_vec_scalar_mul_nmod(Q, A, lenA, invL, mod);
    }
    else
    {
        Q[0] = n_mulmod2_preinv(A[lenA-1], invL, mod.n, mod.ninv);

        _nmod_vec_scalar_mul_nmod(R, B, lenA - 1, Q[0], mod);
        _nmod_vec_sub(R, A, R, lenA - 1, mod);
    }
}

