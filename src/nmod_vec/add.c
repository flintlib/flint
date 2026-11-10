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

void _nmod_vec_add(nn_ptr res, nn_srcptr vec1,
                   nn_srcptr vec2, slong len, nmod_t mod)
{
    slong i;

#if defined(__clang__)
    /* In place, clang checks at run time that the output does not overlap
       the inputs; when res equals an input that check fails and the loop
       runs in scalar code (Nov 2026: 2.5-5.5x slower with clang 18, 19,
       21, 22 and Apple clang 17). gcc vectorizes the in-place case
       already: there the extra loops would only add a cost (up to 2 ns per
       call on short vectors). */
    if (res == vec2)
        FLINT_SWAP(nn_srcptr, vec1, vec2);
#endif

    if (mod.norm)
    {
#if defined(__clang__)
        if (res == vec1)
            for (i = 0 ; i < len; i++)
                res[i] = _nmod_add(res[i], vec2[i], mod);
        else
#endif
        for (i = 0 ; i < len; i++)
            res[i] = _nmod_add(vec1[i], vec2[i], mod);
    } else
    {
#if defined(__clang__)
        if (res == vec1)
            for (i = 0 ; i < len; i++)
                res[i] = nmod_add(res[i], vec2[i], mod);
        else
#endif
        for (i = 0 ; i < len; i++)
            res[i] = nmod_add(vec1[i], vec2[i], mod);
    }
}
