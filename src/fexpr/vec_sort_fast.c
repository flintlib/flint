/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "fexpr.h"

void
_fexpr_vec_sort_fast(fexpr_ptr vec, slong len)
{
    qsort(vec, len, sizeof(fexpr_struct), (int(*)(const void*,const void*)) fexpr_cmp_fast);
}
