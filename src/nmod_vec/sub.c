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

void _nmod_vec_sub(nn_ptr res, nn_srcptr vec1,
                   nn_srcptr vec2, slong len, nmod_t mod)
{
    slong i;

    /* in place with clang, see _nmod_vec_add */
    if (mod.norm)
    {
#if defined(__clang__)
        if (res == vec1)
            for (i = 0 ; i < len; i++)
                res[i] = _nmod_sub(res[i], vec2[i], mod);
        else if (res == vec2)
            for (i = 0 ; i < len; i++)
                res[i] = _nmod_sub(vec1[i], res[i], mod);
        else
#endif
        for (i = 0 ; i < len; i++)
            res[i] = _nmod_sub(vec1[i], vec2[i], mod);
    } else
    {
#if defined(__clang__)
        if (res == vec1)
            for (i = 0 ; i < len; i++)
                res[i] = nmod_sub(res[i], vec2[i], mod);
        else if (res == vec2)
            for (i = 0 ; i < len; i++)
                res[i] = nmod_sub(vec1[i], res[i], mod);
        else
#endif
        for (i = 0 ; i < len; i++)
            res[i] = nmod_sub(vec1[i], vec2[i], mod);
    }
}
