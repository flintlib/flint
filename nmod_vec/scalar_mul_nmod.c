/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2021 Fredrik Johansson

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

void _nmod_vec_scalar_mul_nmod_fullword(mp_ptr res, mp_srcptr vec, 
                               slong len, mp_limb_t c, nmod_t mod)
{
    slong i;

    for (i = 0; i < len; i++)
        NMOD_MUL_FULLWORD(res[i], vec[i], c, mod);
}

void _nmod_vec_scalar_mul_nmod_generic(mp_ptr res, mp_srcptr vec, 
                               slong len, mp_limb_t c, nmod_t mod)
{
    slong i;

    for (i = 0; i < len; i++)
        NMOD_MUL_PRENORM(res[i], vec[i], c << mod.norm, mod);
}

void _nmod_vec_scalar_mul_nmod(mp_ptr res, mp_srcptr vec, 
                               slong len, mp_limb_t c, nmod_t mod)
{
    if (NMOD_BITS(mod) == FLINT_BITS)
        _nmod_vec_scalar_mul_nmod_fullword(res, vec, len, c, mod);
    else if (len > 10)
        _nmod_vec_scalar_mul_nmod_shoup(res, vec, len, c, mod);
    else 
        _nmod_vec_scalar_mul_nmod_generic(res, vec, len, c, mod);
}
