/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

#include "ulong_extras.h"

void TEMPLATE(T, poly_add_si)(
    TEMPLATE(T, poly_t a),
    const TEMPLATE(T, poly_t) b,
    slong c,
    const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, poly_t) t;
    TEMPLATE(T, poly_init)(t, ctx);
    TEMPLATE(T, poly_fit_length)(t, 1, ctx);
    TEMPLATE(T, set_si)(t->coeffs + 0, c, ctx);
    t->length = !TEMPLATE(T, is_zero)(t->coeffs + 0, ctx);
    TEMPLATE(T, poly_add)(a, b, t, ctx);
    TEMPLATE(T, poly_clear)(t, ctx);
}

#endif
