/*
    Copyright (C) 2016 Vincent Delecroix

    This file is part of FLINT

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_vec.h"

typedef int(*__compar_fn_t) (const void *, const void *);
void _fmpq_vec_sort(fmpq * vec, slong len)
{
    qsort(vec, len, sizeof(fmpq), (__compar_fn_t) fmpq_cmp);
}
