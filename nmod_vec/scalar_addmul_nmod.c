/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "nmod.h"

void _nmod_vec_scalar_addmul_nmod_fullword(ulong_ptr res, ulong_srcptr vec, 
				             slong len, ulong c, nmod_t mod)
{
    slong i;
    ulong t;

    for (i = 0; i < len; i++)
    {
        NMOD_MUL_FULLWORD(t, vec[i], c, mod);
        res[i] = nmod_add(res[i], t, mod);
    }
}

void _nmod_vec_scalar_addmul_nmod_generic(ulong_ptr res, ulong_srcptr vec, 
				             slong len, ulong c, nmod_t mod)
{
    slong i;
    ulong t;

    for (i = 0; i < len; i++)
    {
        NMOD_MUL_PRENORM(t, vec[i], c << mod.norm, mod);
        res[i] = _nmod_add(res[i], t, mod);
    }
}

void _nmod_vec_scalar_addmul_nmod_shoup(ulong_ptr res, ulong_srcptr vec, 
				             slong len, ulong c, nmod_t mod)
{
    slong i;
    ulong t, cinv;

    cinv = n_mulmod_precomp_shoup(c, mod.n);

    for (i = 0; i < len; i++)
    {
        t = n_mulmod_shoup(c, vec[i], cinv, mod.n);
        res[i] = _nmod_add(res[i], t, mod);
    }
}

void _nmod_vec_scalar_addmul_nmod(ulong_ptr res, ulong_srcptr vec, 
				             slong len, ulong c, nmod_t mod)
{
    if (NMOD_BITS(mod) == FLINT_BITS)
        _nmod_vec_scalar_addmul_nmod_fullword(res, vec, len, c, mod);
    else if (len > 10)
        _nmod_vec_scalar_addmul_nmod_shoup(res, vec, len, c, mod);
    else
        _nmod_vec_scalar_addmul_nmod_generic(res, vec, len, c, mod);
}
