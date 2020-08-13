/*
    Copyright (C) 2019 Tommy Hofmann

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"
#include "nmod_poly.h"

void
fmpq_poly_set_nmod_poly(fmpq_poly_t rop, const nmod_poly_t op)
{
    slong len = op->length;

    if (len == 0)
    {
        fmpq_poly_zero(rop);
    }
    else
    {
        slong i;
        fmpz_one(rop->den);
        fmpq_poly_fit_length(rop, len);
        for (i = 0; i < len; i++)
            fmpz_set_ui_smod(rop->coeffs + i, op->coeffs[i], op->mod.n);
        _fmpq_poly_set_length(rop, len);
    }
}
