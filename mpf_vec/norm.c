/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpf_vec.h"

void
_mpf_vec_norm(mpf_t res, const mpf * vec, slong len)
{
    slong i;
    mpf_t tmp;
    mpf_init(tmp);

    flint_mpf_set_ui(res, 0);
    for (i = 0; i < len; i++)
    {
        mpf_mul(tmp, vec + i, vec + i);
        mpf_add(res, res, tmp);
    }

    mpf_clear(tmp);
}
