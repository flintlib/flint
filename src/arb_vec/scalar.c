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

void
_arb_vec_scalar_mul(arb_ptr res, arb_srcptr vec, slong len, const arb_t c, slong prec)
{
    slong i;
    for (i = 0; i < len; i++)
        arb_mul(res + i, vec + i, c, prec);
}

void
_arb_vec_scalar_div(arb_ptr res, arb_srcptr vec, slong len, const arb_t c, slong prec)
{
    slong i;
    for (i = 0; i < len; i++)
        arb_div(res + i, vec + i, c, prec);
}

void
_arb_vec_scalar_mul_fmpz(arb_ptr res, arb_srcptr vec, slong len, const fmpz_t c, slong prec)
{
    slong i;
    arf_t t;
    arf_init(t);
    arf_set_fmpz(t, c);
    for (i = 0; i < len; i++)
        arb_mul_arf(res + i, vec + i, t, prec);
    arf_clear(t);
}

void
_arb_vec_scalar_mul_2exp_si(arb_ptr res, arb_srcptr src, slong len, slong c)
{
    slong i;
    for (i = 0; i < len; i++)
        arb_mul_2exp_si(res + i, src + i, c);
}

void
_arb_vec_scalar_addmul(arb_ptr res, arb_srcptr vec, slong len, const arb_t c, slong prec)
{
    slong i;
    for (i = 0; i < len; i++)
        arb_addmul(res + i, vec + i, c, prec);
}
