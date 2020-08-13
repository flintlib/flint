/*
    Copyright (C) 2010 William Hart

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

void _nmod_vec_neg(mp_ptr res, mp_srcptr vec, slong len, nmod_t mod)
{
    slong i;
    for (i = 0 ; i < len; i++)
        res[i] = nmod_neg(vec[i], mod);
}
