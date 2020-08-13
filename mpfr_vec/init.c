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

flint_mpfr *
_mpfr_vec_init(slong length, flint_bitcnt_t prec)
{
    slong i;

    __mpfr_struct *vec =
        (__mpfr_struct *) flint_malloc(length * sizeof(__mpfr_struct));

    for (i = 0; i < length; i++)
        mpfr_init2(vec + i, prec);

    return vec;
}
