/*
    Copyright (C) 2011 Andy Novocin
    Copyright (C) 2011 Sebastian Pancratz

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

void fmpz_poly_hensel_lift_tree(slong *link, fmpz_poly_t *v, fmpz_poly_t *w, 
    fmpz_poly_t f, slong r, const fmpz_t p, slong e0, slong e1, slong inv)
{
    fmpz_t p0, p1;

    fmpz_init(p0);
    fmpz_init(p1);

    fmpz_pow_ui(p0, p, e0);
    fmpz_pow_ui(p1, p, e1 - e0);

    fmpz_poly_hensel_lift_tree_recursive(link, v, w, f, 2*r - 4, inv, p0, p1);

    fmpz_clear(p0);
    fmpz_clear(p1);
}

