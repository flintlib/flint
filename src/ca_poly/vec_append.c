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
ca_poly_vec_append(ca_poly_vec_t vec, const ca_poly_t f, ca_ctx_t ctx)
{
    _ca_poly_vec_fit_length(vec, vec->length + 1, ctx);
    ca_poly_set(vec->entries + vec->length, f, ctx);
    vec->length++;
}
