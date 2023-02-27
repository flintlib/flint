/*
    Copyright (C) 2017 Luca De Feo

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

FLINT_DLL void TEMPLATE(T, poly_set_fmpz_mod_poly)(TEMPLATE(T, poly_t) rop,
                                                   const fmpz_mod_poly_t op,
                                                   const TEMPLATE(T, ctx_t) ctx)
{
    slong i, len = op->length;

    TEMPLATE(T, poly_fit_length)(rop, len, ctx);
    _TEMPLATE(T, poly_set_length)(rop, len, ctx);

    for (i = 0; i < len; i++)
        TEMPLATE(T, set_fmpz)(rop->coeffs + i, op->coeffs + i, ctx);
}


#endif
