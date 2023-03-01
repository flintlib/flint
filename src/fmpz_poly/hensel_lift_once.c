/*
    Copyright (C) 2011 Andy Novocin
    Copyright (C) 2011 William Hart
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "fmpz_poly.h"
#include "nmod_poly.h"

void fmpz_poly_hensel_lift_once(fmpz_poly_factor_t lifted_fac, 
                                const fmpz_poly_t f, 
                                const nmod_poly_factor_t local_fac, slong N)
{
    const slong r = local_fac->num;

    slong i;
    slong *link;
    fmpz_poly_t *v, *w;

    link = flint_malloc((2*r - 2) * sizeof(slong));
    v    = flint_malloc(2*(2*r - 2) * sizeof(fmpz_poly_t));
    w    = v + (2*r - 2);

    for(i = 0; i < 2*r - 2; i++)
    {
        fmpz_poly_init(v[i]);
        fmpz_poly_init(w[i]);
    }

    _fmpz_poly_hensel_start_lift(lifted_fac, link, v, w, f, local_fac, N);

    for (i = 0; i < 2*r - 2; i++)
    {
        fmpz_poly_clear(v[i]);
        fmpz_poly_clear(w[i]);
    }
    flint_free(link);
    flint_free(v);
}

