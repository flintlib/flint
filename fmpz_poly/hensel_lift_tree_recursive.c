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

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

void fmpz_poly_hensel_lift_tree_recursive(len_t *link, 
    fmpz_poly_t *v, fmpz_poly_t *w, fmpz_poly_t f, len_t j, len_t inv, 
    const fmpz_t p0, const fmpz_t p1)
{
    if (j >= 0)
    {
        if (inv == 1)
            fmpz_poly_hensel_lift(v[j], v[j + 1], w[j], w[j + 1], f, 
                                  v[j], v[j + 1], w[j], w[j + 1], 
                                  p0, p1);
        else if (inv == -1)
            fmpz_poly_hensel_lift_only_inverse(w[j], w[j+1], 
                                 v[j], v[j+1], w[j], w[j+1], p0, p1);
        else
            fmpz_poly_hensel_lift_without_inverse(v[j], v[j+1], f, 
                                                  v[j], v[j+1], w[j], w[j+1], 
                                                  p0, p1);

        fmpz_poly_hensel_lift_tree_recursive(link, v, w, v[j], link[j], 
            inv, p0, p1);
        fmpz_poly_hensel_lift_tree_recursive(link, v, w, v[j+1], link[j+1], 
            inv, p0, p1);
    }
}

