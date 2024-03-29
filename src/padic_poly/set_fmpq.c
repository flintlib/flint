/*
    Copyright (C) 2011, 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "padic.h"
#include "padic_poly.h"

void padic_poly_set_fmpq(padic_poly_t poly, const fmpq_t x,
                         const padic_ctx_t ctx)
{
    padic_t y;

    padic_init2(y, padic_poly_prec(poly));
    padic_set_fmpq(y, x, ctx);
    padic_poly_set_padic(poly, y, ctx);
    padic_clear(y);
}
