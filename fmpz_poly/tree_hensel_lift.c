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
    Computes $p_0 = p^{e_0}$ and $p_1 = p^{e_1 - e_0}$ for a small prime $p$ 
    and $P = p^{e_1}$.

    If we aim to lift to $p^b$ then $f$ is the polynomial whose factors we 
    wish to lift, made monic mod $p^b$. As usual, \code{{link, v, w}} is an 
    initialised tree.

    This starts the recursion on lifting the \emph{product tree} for lifting 
    from $p^{e_0} to $p^{e_1}$. The value of \code{inv} corresponds to that 
    given for the function \code{fmpz_poly_tree_hensel_lift_recursive()}. We 
    set $r$ to the number of local factors of $f$.

    In terms of the notation, above $P = p^{e_1}$, $p_0 = p^{e_0}$ and 
    $p_1 = p^{e_1-e_0}$.

    Assumes that $f$ is monic.
 */

void fmpz_poly_tree_hensel_lift(long *link, fmpz_poly_t *v, fmpz_poly_t *w, fmpz_t P, 
    fmpz_poly_t f, long r, const fmpz_t p, long e0, long e1, long inv)
{
    fmpz_t p0, p1;

    fmpz_init(p0);
    fmpz_init(p1);

    fmpz_pow_ui(p0, p, e0);
    fmpz_pow_ui(p1, p, e1 - e0);

    fmpz_mul(P, p0, p1);

    fmpz_poly_tree_hensel_lift_recursive(link, v, w, f, 2*r - 4, inv, p0, p1, P);

    fmpz_clear(p0);
    fmpz_clear(p1);
}

