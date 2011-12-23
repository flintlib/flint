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

void fmpz_poly_hensel_lift_only_inverse(fmpz_poly_t Aout, fmpz_poly_t Bout, 
    const fmpz_poly_t f, 
    const fmpz_poly_t G, const fmpz_poly_t H, 
    const fmpz_poly_t a, const fmpz_poly_t b, 
    const fmpz_t p, const fmpz_t p1)
{
    fmpz_poly_t A, B;
    fmpz_poly_t a1, b1, t1, t2, r, unity;
    fmpz_mod_poly_t gg, hh, aa, bb, g, h;
    fmpz_mod_poly_t rr, aa1, bb1;

    fmpz_poly_init(A);
    fmpz_poly_init(B);

    fmpz_poly_init(a1);
    fmpz_poly_init(b1);
    fmpz_poly_init(t1);
    fmpz_poly_init(t2);
    fmpz_poly_init(r);
    fmpz_poly_init(unity);

    /* Make a check that c is divisible by p */
    fmpz_mod_poly_init(gg, p1);
    fmpz_mod_poly_init(hh, p1);
    fmpz_mod_poly_init(g, p);
    fmpz_mod_poly_init(h, p);
    fmpz_mod_poly_init(aa, p1);
    fmpz_mod_poly_init(bb, p1);

    fmpz_mod_poly_init(rr, p1);
    fmpz_mod_poly_init(aa1, p1);
    fmpz_mod_poly_init(bb1, p1);

    fmpz_mod_poly_set_fmpz_poly(g, G); /*reduce mod p*/
    fmpz_mod_poly_set_fmpz_poly(h, H); /*reduce mod p*/
    fmpz_mod_poly_set(gg, g); /*embed mod p1*/
    fmpz_mod_poly_set(hh, h); /*embed mod p1*/
    fmpz_mod_poly_set_fmpz_poly(aa, a);
    fmpz_mod_poly_set_fmpz_poly(bb, b);

    /*Lifting the inverses now*/
    fmpz_poly_set_coeff_si(unity, 0, -1);
    fmpz_poly_mul(t1, a, G);
    fmpz_poly_mul(t2, b, H);

    fmpz_poly_add(t1, t1, t2);
    fmpz_poly_add(t1, t1, unity);
    fmpz_poly_neg(t1, t1);

    fmpz_poly_scalar_divexact_fmpz(r, t1, p);
    fmpz_mod_poly_set_fmpz_poly(rr, r);

    fmpz_mod_poly_rem(bb1, rr, gg);
    fmpz_mod_poly_mul(bb1, bb1, bb);
    fmpz_mod_poly_rem(bb1, bb1, gg);

    fmpz_mod_poly_rem(aa1, rr, hh);
    fmpz_mod_poly_mul(aa1, aa1, aa);
    fmpz_mod_poly_rem(aa1, aa1, hh);

    fmpz_mod_poly_get_fmpz_poly(a1, aa1);
    fmpz_poly_scalar_mul_fmpz(a1, a1, p);
    fmpz_poly_add(A, a, a1);

    fmpz_mod_poly_get_fmpz_poly(b1, bb1);
    fmpz_poly_scalar_mul_fmpz(b1, b1, p);
    fmpz_poly_add(B, b, b1);

    fmpz_poly_set(Aout, A);
    fmpz_poly_set(Bout, B);

    fmpz_mod_poly_clear(rr);
    fmpz_mod_poly_clear(aa1);
    fmpz_mod_poly_clear(bb1);

    fmpz_poly_clear(a1);
    fmpz_poly_clear(b1);
    fmpz_poly_clear(t1);
    fmpz_poly_clear(t2);
    fmpz_poly_clear(r);
    fmpz_poly_clear(unity);

    fmpz_mod_poly_clear(gg);
    fmpz_mod_poly_clear(hh);
    fmpz_mod_poly_clear(aa);
    fmpz_mod_poly_clear(bb);
    fmpz_mod_poly_clear(g);
    fmpz_mod_poly_clear(h);

    fmpz_poly_clear(A);
    fmpz_poly_clear(B);
}

