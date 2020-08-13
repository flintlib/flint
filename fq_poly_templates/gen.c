/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

void
TEMPLATE(T, poly_gen) (TEMPLATE(T, poly_t) f, const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, poly_fit_length) (f, 2, ctx);
    TEMPLATE(T, zero) (f->coeffs, ctx);
    TEMPLATE(T, one) (f->coeffs + 1, ctx);
    _TEMPLATE(T, poly_set_length) (f, 2, ctx);
}


#endif
