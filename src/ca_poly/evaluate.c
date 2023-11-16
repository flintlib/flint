/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_poly.h"

void
_ca_poly_evaluate(ca_t res, ca_srcptr f, slong len,
                           const ca_t x, ca_ctx_t ctx)
{
    _ca_poly_evaluate_horner(res, f, len, x, ctx);
}

void
ca_poly_evaluate(ca_t res, const ca_poly_t f, const ca_t a, ca_ctx_t ctx)
{
    _ca_poly_evaluate(res, f->coeffs, f->length, a, ctx);
}
