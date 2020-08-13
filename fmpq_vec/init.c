/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2016 Vincent Delecroix

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "fmpq.h"
#include "fmpq_vec.h"

fmpq *
_fmpq_vec_init(slong len)
{
    fmpq * v = (fmpq *) flint_malloc(sizeof(fmpq) * len);
    slong i;

    for (i = 0; i < len; i++)
        fmpq_init(v + i);

    return v;
}
