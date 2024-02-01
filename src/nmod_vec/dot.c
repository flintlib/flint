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

mp_limb_t
_nmod_vec_dot(mp_srcptr vec1, mp_srcptr vec2, slong len, nmod_t mod, int nlimbs)
{
    mp_limb_t res;
    slong i;
    NMOD_VEC_DOT(res, i, len, vec1[i], vec2[i], mod, nlimbs);
    return res;
}

int
_nmod_vec_dot_bound_limbs(slong len, nmod_t mod)
{
    mp_limb_t t2, t1, t0, u1, u0;

    umul_ppmm(t1, t0, mod.n - 1, mod.n - 1);
    umul_ppmm(t2, t1, t1, len);
    umul_ppmm(u1, u0, t0, len);
    add_sssaaaaaa(t2, t1, t0,  t2, t1, UWORD(0),  UWORD(0), u1, u0);

    if (t2 != 0) return 3;
    if (t1 != 0) return 2;
    return (t0 != 0);
}

mp_limb_t
_nmod_vec_dot_ptr(mp_srcptr vec1, const mp_ptr * vec2, slong offset,
                            slong len, nmod_t mod, int nlimbs)
{
    mp_limb_t res;
    slong i;
    NMOD_VEC_DOT(res, i, len, vec1[i], vec2[i][offset], mod, nlimbs);
    return res;
}

static mp_limb_t
nmod_fmma(mp_limb_t a, mp_limb_t b, mp_limb_t c, mp_limb_t d, nmod_t mod)
{
    a = nmod_mul(a, b, mod);
    NMOD_ADDMUL(a, c, d, mod);
    return a;
}

mp_limb_t
_nmod_vec_dot_rev(mp_srcptr vec1, mp_srcptr vec2, slong len, nmod_t mod, int nlimbs)
{
    mp_limb_t res;
    slong i;

    if (len <= 2 && nlimbs >= 2)
    {
        if (len == 2)
            return nmod_fmma(vec1[0], vec2[1], vec1[1], vec2[0], mod);
        if (len == 1)
            return nmod_mul(vec1[0], vec2[0], mod);
        return 0;
    }

    NMOD_VEC_DOT(res, i, len, vec1[i], vec2[len - 1 - i], mod, nlimbs);
    return res;
}
