/*
    Copyright (C) 2015 Tommy Hofmann

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <stdio.h>
#undef ulong
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fq.h"
#include "fq_poly.h"

void fq_poly_factor_get_poly(fq_poly_t z,
                       const fq_poly_factor_t fac, slong i, const fq_ctx_t ctx)
{
    fq_poly_set(z, fac->poly + i, ctx);
}

