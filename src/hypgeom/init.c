/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"
#include "mag.h"
#include "hypgeom.h"

void
hypgeom_init(hypgeom_t hyp)
{
    fmpz_poly_init(hyp->A);
    fmpz_poly_init(hyp->B);
    fmpz_poly_init(hyp->P);
    fmpz_poly_init(hyp->Q);
    mag_init(hyp->MK);
    hyp->have_precomputed = 0;
}

void
hypgeom_clear(hypgeom_t hyp)
{
    fmpz_poly_clear(hyp->A);
    fmpz_poly_clear(hyp->B);
    fmpz_poly_clear(hyp->P);
    fmpz_poly_clear(hyp->Q);
    mag_clear(hyp->MK);
}
