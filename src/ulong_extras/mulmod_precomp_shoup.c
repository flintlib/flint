/*
    Copyright (C) 2015 William Hart
    Copyright (C) 2015 Vladimir Glazachev

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "ulong_extras.h"

/* Computes the W' = [w * b / p] (b = ulong power) */
ulong
n_mulmod_precomp_shoup(ulong w, ulong p)
{
   ulong q, r;
   udiv_qrnnd(q, r, w, UWORD(0), p);
   return q;
}
