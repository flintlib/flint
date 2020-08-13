/*
    Copyright (C) 2015 William Hart
    Copyright (C) 2015 Vladimir Glazachev

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include <stdlib.h>
#include "flint.h"
#include "ulong_extras.h"
#include "nmod_vec.h"

void _nmod_vec_scalar_mul_nmod_shoup(mp_ptr res, mp_srcptr vec, 
                               slong len, mp_limb_t c, nmod_t mod)
{
    slong i;
    mp_limb_t w_pr;
    w_pr = n_mulmod_precomp_shoup(c, mod.n);
    for (i = 0; i < len; i++)
        res[i] = n_mulmod_shoup(c, vec[i], w_pr, mod.n);
}
