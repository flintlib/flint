/*
    Copyright (C) 2011 Andy Novocin
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

void fmpz_poly_hensel_lift_tree_recursive(slong *link,
    fmpz_poly_t *v, fmpz_poly_t *w, fmpz_poly_t f, slong j, slong inv,
    const fmpz_t p0, const fmpz_t p1)
{
    flint_printf("enter fmpz_poly_hensel_lift_tree_recursive\n"); fflush(stdout);
    if (j >= 0)
    {
        if (inv == 1)
        {
            flint_printf("mid1 fmpz_poly_hensel_lift_tree_recursive\n"); fflush(stdout);
            fmpz_poly_hensel_lift(v[j], v[j + 1], w[j], w[j + 1], f,
                                  v[j], v[j + 1], w[j], w[j + 1],
                                  p0, p1);
        }
        else if (inv == -1)
        {
            flint_printf("mid2 fmpz_poly_hensel_lift_tree_recursive\n"); fflush(stdout);
            fmpz_poly_hensel_lift_only_inverse(w[j], w[j+1],
                                 v[j], v[j+1], w[j], w[j+1], p0, p1);
        }
        else
        {
            flint_printf("mid3 fmpz_poly_hensel_lift_tree_recursive\n"); fflush(stdout);
            fmpz_poly_hensel_lift_without_inverse(v[j], v[j+1], f,
                                                  v[j], v[j+1], w[j], w[j+1],
                                                  p0, p1);
        }

        flint_printf("mid4 fmpz_poly_hensel_lift_tree_recursive\n"); fflush(stdout);
        fmpz_poly_hensel_lift_tree_recursive(link, v, w, v[j], link[j],
            inv, p0, p1);
        flint_printf("mid5 fmpz_poly_hensel_lift_tree_recursive\n"); fflush(stdout);
        fmpz_poly_hensel_lift_tree_recursive(link, v, w, v[j+1], link[j+1],
            inv, p0, p1);
    }
    flint_printf("exit fmpz_poly_hensel_lift_tree_recursive\n"); fflush(stdout);
}

