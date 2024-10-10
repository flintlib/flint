/*
    Copyright (C) 2010, 2015 William Hart
    Copyright (C) 2015 Vladimir Glazachev
    Copyright (C) 2021 Fredrik Johansson
    Copyright (C) 2024 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint-mparam.h"
#include "nmod_vec.h"

static void _nmod_vec_scalar_addmul_nmod_fullword(nn_ptr res, nn_srcptr vec,
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

void _nmod_vec_scalar_addmul_nmod_generic(nn_ptr res, nn_srcptr vec,
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

void _nmod_vec_scalar_addmul_nmod_shoup(nn_ptr res, nn_srcptr vec,
				             slong len, ulong c, nmod_t mod)
{
    const ulong c_pr = n_mulmod_precomp_shoup(c, mod.n);

    for (slong i = 0; i < len; i++)
    {
        ulong t = n_mulmod_shoup(c, vec[i], c_pr, mod.n);
        res[i] = _nmod_add(res[i], t, mod);
    }
}

void _nmod_vec_scalar_addmul_nmod(nn_ptr res, nn_srcptr vec,
				             slong len, ulong c, nmod_t mod)
{
    if (c == UWORD(0))
        return;
    else if (c == UWORD(1))
        _nmod_vec_add(res, res, vec, len, mod);
    else if (c == mod.n - UWORD(1))
        _nmod_vec_sub(res, res, vec, len, mod);
    else if (NMOD_BITS(mod) == FLINT_BITS)
        _nmod_vec_scalar_addmul_nmod_fullword(res, vec, len, c, mod);
    else if (len >= FLINT_MULMOD_SHOUP_THRESHOLD)
        _nmod_vec_scalar_addmul_nmod_shoup(res, vec, len, c, mod);
    else
        _nmod_vec_scalar_addmul_nmod_generic(res, vec, len, c, mod);
}

static void _nmod_vec_scalar_mul_nmod_fullword(nn_ptr res, nn_srcptr vec,
                               slong len, ulong c, nmod_t mod)
{
    slong i;

    for (i = 0; i < len; i++)
        NMOD_MUL_FULLWORD(res[i], vec[i], c, mod);
}

void _nmod_vec_scalar_mul_nmod_generic(nn_ptr res, nn_srcptr vec,
                               slong len, ulong c, nmod_t mod)
{
    slong i;

    for (i = 0; i < len; i++)
        NMOD_MUL_PRENORM(res[i], vec[i], c << mod.norm, mod);
}

void _nmod_vec_scalar_mul_nmod_shoup(nn_ptr res, nn_srcptr vec,
                               slong len, ulong c, nmod_t mod)
{
    const ulong c_pr = n_mulmod_precomp_shoup(c, mod.n);

    for (slong i = 0; i < len; i++)
        res[i] = n_mulmod_shoup(c, vec[i], c_pr, mod.n);
}

void _nmod_vec_scalar_mul_nmod(nn_ptr res, nn_srcptr vec,
                               slong len, ulong c, nmod_t mod)
{
    if (c == UWORD(0))
        _nmod_vec_zero(res, len);
    else if (c == UWORD(1))
        _nmod_vec_set(res, vec, len);
    else if (c == mod.n - UWORD(1))
        _nmod_vec_neg(res, vec, len, mod);
    else if (NMOD_BITS(mod) == FLINT_BITS)
        _nmod_vec_scalar_mul_nmod_fullword(res, vec, len, c, mod);
    else if (len >= FLINT_MULMOD_SHOUP_THRESHOLD)
        _nmod_vec_scalar_mul_nmod_shoup(res, vec, len, c, mod);
    else
        _nmod_vec_scalar_mul_nmod_generic(res, vec, len, c, mod);
}
