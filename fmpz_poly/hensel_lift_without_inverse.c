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
#include "fmpz_mod_poly.h"

#define lift(G, g, b)                                                                     \
do {                                                                                      \
    fmpz_poly_fit_length(G, g->length);                                                   \
    _fmpz_poly_set_length(G, g->length);                                                  \
    _fmpz_mod_poly_mul(T1, C, f->length, b->coeffs, b->length, p1);                       \
    _fmpz_mod_poly_rem(T2, T1, f->length + b->length - 1, g->coeffs, g->length, one, p1); \
    _fmpz_vec_scalar_mul_fmpz(T1, T2, g->length, p);                                      \
    _fmpz_vec_add(G->coeffs, T1, g->coeffs, g->length);                                   \
    _fmpz_poly_normalise(G);                                                              \
} while(0);


/*
    This is the main Hensel lifting routine, which performs a Hensel step
    from polynomials mod $p$ to polynomials mod $p p_1$.

    One starts with polynomials $f$, $g$, $h$ such that $f = gh \pmod{p}$. 
    The polynomials $a$, $b$ satisfy $ag + bh = 1 \pmod{p}$.

    Upon return we have $A G + B H = 1 \pmod{p p_1}$ and 
    $f = G H \pmod{p p_1}$, where $g = G \pmod{p}$, etc.

    We require that $1 < p_1 \leq p$.

    The output arguments $G$ and $H$ may only be aliased with 
    the input arguments $g$ and $h$, respectively.
 */

void fmpz_poly_hensel_lift_without_inverse(fmpz_poly_t G, fmpz_poly_t H, 
    const fmpz_poly_t f, const fmpz_poly_t g, const fmpz_poly_t h, 
    const fmpz_poly_t a, const fmpz_poly_t b, 
    const fmpz_t p, const fmpz_t p1)
{
    fmpz one[1] = {1L};
    fmpz *C, *T1, *T2;
    const long lenT = f->length + FLINT_MAX(a->length, b->length) - 1;
    const long len  = f->length + 2 * lenT;

    C  = _fmpz_vec_init(len);
    T1 = C  + f->length;
    T2 = T1 + lenT;

    /* {C, lenF} := {(f-gh)/p, lenF} */
    if (g->length >= h->length)
        _fmpz_poly_mul(C, g->coeffs, g->length, h->coeffs, h->length);
    else
        _fmpz_poly_mul(C, h->coeffs, h->length, g->coeffs, g->length);
    _fmpz_poly_sub(C, f->coeffs, f->length, C, f->length);
    _fmpz_vec_scalar_divexact_fmpz(C, C, f->length, p);

    lift(G, g, b);
    lift(H, h, a);

    _fmpz_vec_clear(C, len);
}

