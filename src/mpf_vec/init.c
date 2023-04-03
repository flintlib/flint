/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpf_vec.h"

mpf *
_mpf_vec_init(slong len, flint_bitcnt_t prec)
{
    slong i;

    mpf *vec = (mpf *) flint_malloc(len * sizeof(mpf));

    for (i = 0; i < len; i++)
        mpf_init2(vec + i, prec);

    return vec;
}
