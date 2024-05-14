/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "arb_vec.h"

arb_ptr _arb_vec_init(slong n)
{
    slong i;
    arb_ptr v = flint_malloc(sizeof(arb_struct) * n);

    for (i = 0; i < n; i++)
        arb_init(v + i);

    return v;
}

void _arb_vec_clear(arb_ptr v, slong n)
{
    slong i;
    for (i = 0; i < n; i++)
        arb_clear(v + i);
    flint_free(v);
}

arb_ptr _arb_vec_entry_ptr(arb_ptr vec, slong i)
{
    return vec + i;
}

void
_arb_vec_swap(arb_ptr res, arb_ptr vec, slong len)
{
    slong i;
    for (i = 0; i < len; i++)
        arb_swap(res + i, vec + i);
}

slong _arb_vec_allocated_bytes(arb_srcptr vec, slong len)
{
    slong i, size;

    size = len * sizeof(arb_struct);

    for (i = 0; i < len; i++)
        size += arb_allocated_bytes(vec + i);

    return size;
}

double _arb_vec_estimate_allocated_bytes(slong len, slong prec)
{
    double size;

    size = len * (double) sizeof(arb_struct);

    if (prec > ARF_NOPTR_LIMBS * FLINT_BITS)
        size += len * (double) ((prec + FLINT_BITS - 1) / FLINT_BITS) * sizeof(ulong);

    return size;
}

void _arb_vec_trim(arb_ptr res, arb_srcptr vec, slong len)
{
    slong i;
    for (i = 0; i < len; i++)
        arb_trim(res + i, vec + i);
}

slong _arb_vec_bits(arb_srcptr x, slong len)
{
    slong i, b, c;

    b = 0;
    for (i = 0; i < len; i++)
    {
        c = arb_bits(x + i);
        b = FLINT_MAX(b, c);
    }

    return b;
}
