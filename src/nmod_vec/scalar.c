/*
    Copyright (C) 2010, 2015 William Hart
    Copyright (C) 2015 Vladimir Glazachev
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "nmod_vec.h"

static void _nmod_vec_scalar_addmul_nmod_fullword(mp_ptr res, mp_srcptr vec,
				             slong len, mp_limb_t c, nmod_t mod)
{
    slong i;
    mp_limb_t t;

    for (i = 0; i < len; i++)
    {
        NMOD_MUL_FULLWORD(t, vec[i], c, mod);
        res[i] = nmod_add(res[i], t, mod);
    }
}

static void _nmod_vec_scalar_addmul_nmod_generic(mp_ptr res, mp_srcptr vec,
				             slong len, mp_limb_t c, nmod_t mod)
{
    slong i;
    mp_limb_t t;

    for (i = 0; i < len; i++)
    {
        NMOD_MUL_PRENORM(t, vec[i], c << mod.norm, mod);
        res[i] = _nmod_add(res[i], t, mod);
    }
}

static void _nmod_vec_scalar_addmul_nmod_shoup(mp_ptr res, mp_srcptr vec,
				             slong len, mp_limb_t c, nmod_t mod)
{
    slong i;
    mp_limb_t t, cinv;

    cinv = n_mulmod_precomp_shoup(c, mod.n);

    for (i = 0; i < len; i++)
    {
        t = n_mulmod_shoup(c, vec[i], cinv, mod.n);
        res[i] = _nmod_add(res[i], t, mod);
    }
}

void _nmod_vec_scalar_addmul_nmod(mp_ptr res, mp_srcptr vec,
				             slong len, mp_limb_t c, nmod_t mod)
{
    if (NMOD_BITS(mod) == FLINT_BITS)
        _nmod_vec_scalar_addmul_nmod_fullword(res, vec, len, c, mod);
    else if (len > 10)
        _nmod_vec_scalar_addmul_nmod_shoup(res, vec, len, c, mod);
    else
        _nmod_vec_scalar_addmul_nmod_generic(res, vec, len, c, mod);
}

static void _nmod_vec_scalar_mul_nmod_fullword(mp_ptr res, mp_srcptr vec,
                               slong len, mp_limb_t c, nmod_t mod)
{
    slong i;

    for (i = 0; i < len; i++)
        NMOD_MUL_FULLWORD(res[i], vec[i], c, mod);
}

static void _nmod_vec_scalar_mul_nmod_generic(mp_ptr res, mp_srcptr vec,
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

void _nmod_vec_scalar_mul_nmod_shoup(mp_ptr res, mp_srcptr vec,
                               slong len, mp_limb_t c, nmod_t mod)
{
    slong i;
    mp_limb_t w_pr;
    w_pr = n_mulmod_precomp_shoup(c, mod.n);
    for (i = 0; i < len; i++)
        res[i] = n_mulmod_shoup(c, vec[i], w_pr, mod.n);
}
