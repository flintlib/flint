/*
    Copyright (C) 2010 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <gmp.h>
#include <mpfr.h>
#include "flint.h"
#include "mpfr_vec.h"

void
_mpfr_vec_add(mpfr_ptr res, mpfr_srcptr vec1, mpfr_srcptr vec2, slong length)
{
    slong i;
    for (i = 0; i < length; i++)
        mpfr_add(res + i, vec1 + i, vec2 + i, GMP_RNDN);
}
