/*
    Copyright (C) 2015 Vladimir Glazachev

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "aprcl.h"

void
_unity_zpq_mul_unity_p(unity_zpq f)
{
    ulong i;
    for (i = f->p - 1; i > 0; i--)
        fmpz_mod_poly_swap(f->polys[i], f->polys[i - 1]);
}

/*
    Computes unity_zpq * \zeta_p by swapping poly coeffs.
*/
void
unity_zpq_mul_unity_p_pow(unity_zpq f, const unity_zpq g, ulong k)
{
    ulong i;
    unity_zpq_copy(f, g);
    for (i = 0; i < k; i++)
        _unity_zpq_mul_unity_p(f);
}

