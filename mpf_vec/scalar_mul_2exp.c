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

void
_mpf_vec_scalar_mul_2exp(mpf * res, const mpf * vec, slong len,
                         flint_bitcnt_t exp)
{
    slong i;
    for (i = 0; i < len; i++)
        mpf_mul_2exp(res + i, vec + i, exp);
}
