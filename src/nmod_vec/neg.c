/*
    Copyright (C) 2010 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "nmod_vec.h"

void _nmod_vec_neg(nn_ptr res, nn_srcptr vec, slong len, nmod_t mod)
{
    slong i;
    for (i = 0 ; i < len; i++)
        res[i] = nmod_neg(vec[i], mod);
}
