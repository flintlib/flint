/*
    Copyright (C) 2015 Vladimir Glazachev

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "aprcl.h"

slong
unity_zpq_p_unity(const unity_zpq f)
{
    slong i;
    ulong is_punity;

    is_punity = f->p;

    for (i = 0; i < f->p; i++)
    {
        if (fmpz_equal_ui(f->polys[i]->coeffs + i, 1))
        {
            if (is_punity != f->p)
                return 0;

            is_punity = i;
        }
    }

    return is_punity;
}

int
unity_zpq_is_p_unity(const unity_zpq f)
{
    if (unity_zpq_p_unity(f) != f->p)
        return 1;

    return 0;
}

int
unity_zpq_is_p_unity_generator(const unity_zpq f)
{
    slong upow = unity_zpq_p_unity(f);

    if (upow != f->p && upow != 0)
            return 1;

    return 0;
}

