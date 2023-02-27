/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"

#ifndef __compar_fn_t
#if defined(_MSC_VER)
typedef int(*__compar_fn_t) (const void *, const void *);
#else
typedef int(*__compar_fn_t) (__const void *, __const void *);
#endif
#endif

void _fmpz_vec_sort(fmpz * vec, slong len)
{
    qsort(vec, len, sizeof(fmpz), (__compar_fn_t) fmpz_cmp);
}
