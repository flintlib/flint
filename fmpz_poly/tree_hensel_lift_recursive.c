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
    Takes a current Hensel tree {link, v, w} and a pair {j, j+1} 
    of entries in the tree and lifts the tree from mod $p_0$ to 
    mod $P = p_0 p_1$.

    Set \code{inv} to $-1$ if restarting Hensel lifting, $0$ if stopping 
    and $1$ otherwise. 

    Here $f = g h$ is the polynomial whose factors we are trying to lift. 
    We will have that \code{v[j]} is the product of \code{v[link[j]]} and 
    \code{v[link[j] + 1]} as described above.

    Does support aliasing of $f$ with one of the polynomials in 
    the lists $v$ and $w$.  But the polynomials in these two lists 
    are not allowed to be aliases of each other.
 */

void fmpz_poly_tree_hensel_lift_recursive(long *link, 
    fmpz_poly_t *v, fmpz_poly_t *w, fmpz_poly_t f, long j, long inv, 
    const fmpz_t p0, const fmpz_t p1, const fmpz_t P)
{
    if (j >= 0)
    {
        if (inv == 1)
            fmpz_poly_hensel_lift(v[j], v[j + 1], w[j], w[j + 1], f, 
                                  v[j], v[j + 1], w[j], w[j + 1], 
                                  p0, p1);
        else if (inv == -1)
            fmpz_poly_hensel_lift_only_inverse(w[j], w[j+1], f, 
                                 v[j], v[j+1], w[j], w[j+1], p0, p1);
        else
            fmpz_poly_hensel_lift_without_inverse(v[j], v[j+1], f, 
                                                  v[j], v[j+1], w[j], w[j+1], 
                                                  p0, p1);

        fmpz_poly_tree_hensel_lift_recursive(link, v, w, v[j], link[j], 
            inv, p0, p1, P);
        fmpz_poly_tree_hensel_lift_recursive(link, v, w, v[j+1], link[j+1], 
            inv, p0, p1, P);
    }
}

