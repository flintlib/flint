/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2011 Andy Novocin
    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

/*
    XXX: Can alias f and v, w.
 */

void fmpz_poly_tree_hensel_lift_recursive(long *link, 
    fmpz_poly_t *v, fmpz_poly_t *w, fmpz_poly_t f, long j, long inv, 
    const fmpz_t p, const fmpz_t p1, const fmpz_t big_P)
{
    if (j >= 0)
    {
        if (inv == 1)
            fmpz_poly_hensel_lift(v[j], v[j + 1], w[j], w[j + 1], f, 
                                  v[j], v[j + 1], w[j], w[j + 1], 
                                  p, p1, big_P);
        else if (inv == -1)
            fmpz_poly_hensel_lift_only_inverse(w[j], w[j + 1], f, 
                               v[j], v[j + 1], w[j], w[j + 1], p, p1, big_P);
        else
            fmpz_poly_hensel_lift_without_inverse(v[j], v[j + 1], f, 
                                                  v[j], v[j + 1], w[j], w[j + 1], 
                                                  p, p1, big_P);

        fmpz_poly_tree_hensel_lift_recursive(link, v, w, v[j], link[j], inv, p, p1, big_P);
        fmpz_poly_rec_tree_hensel_lift(link, v, w, v[j + 1], link[j + 1], inv, p, p1, big_P);
    }
}

