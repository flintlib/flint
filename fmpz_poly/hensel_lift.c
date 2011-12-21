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

/*
    N.B.  All input polynomials have degree at least 1.

    TODO:  Change in interface, assume that p1 > 1.

f = g h mod p
ag + bh = 1 mod p

G := ([(f - gh)/p] b mod g) p + g
H := ([(f - gh)/p] a mod h) p + h

B := ([(1 - aG - bH)/p] b mod g) p + b
A := ([(1 - aG - bH)/p] a mod h) p + a

 */

void fmpz_poly_hensel_lift(fmpz_poly_t Gout, fmpz_poly_t Hout, 
    fmpz_poly_t Aout, fmpz_poly_t Bout, 
    const fmpz_poly_t f, 
    const fmpz_poly_t g, const fmpz_poly_t h, 
    const fmpz_poly_t a, const fmpz_poly_t b, 
    const fmpz_t p, const fmpz_t p1, const fmpz_t big_P)
{
    fmpz_poly_t c, g1, h1, G, H, A, B;

    fmpz_mod_poly_t cc, gg, hh, aa, bb, tt, gg1, hh1;

    fmpz_poly_t a1, b1, t1, t2, r, unity;
    fmpz_mod_poly_t rr, aa1, bb1;

    fmpz_poly_init(c);
    fmpz_poly_init(g1);
    fmpz_poly_init(h1);
    fmpz_poly_init(G);
    fmpz_poly_init(H);
    fmpz_poly_init(A);
    fmpz_poly_init(B);

    fmpz_poly_mul(c, g, h);
    fmpz_poly_sub(c, f, c);
    /* Make a check that c is divisible by p */
    fmpz_poly_scalar_divexact_fmpz(c, c, p);

        fmpz_mod_poly_init(cc, p1);
        fmpz_mod_poly_init(gg, p1);
        fmpz_mod_poly_init(hh, p1);
        fmpz_mod_poly_init(aa, p1);
        fmpz_mod_poly_init(bb, p1);
        fmpz_mod_poly_init(tt, p1);
        fmpz_mod_poly_init(gg1, p1);
        fmpz_mod_poly_init(hh1, p1);

        fmpz_mod_poly_set_fmpz_poly(cc, c);
        fmpz_mod_poly_set_fmpz_poly(gg, g);
        fmpz_mod_poly_set_fmpz_poly(hh, h);
        fmpz_mod_poly_set_fmpz_poly(aa, a);
        fmpz_mod_poly_set_fmpz_poly(bb, b);

        /*
            TODO: Can save some non-trivial time here by precomputing 
            inverses for gg and hh
            When I make a precomputing function use GG, HH instead
         */

        fmpz_mod_poly_rem(gg1, cc, gg);
        fmpz_mod_poly_mul(gg1, gg1, bb);
        fmpz_mod_poly_rem(gg1, gg1, gg);

        fmpz_mod_poly_rem(hh1, cc, hh);
        fmpz_mod_poly_mul(hh1, hh1, aa);
        fmpz_mod_poly_rem(hh1, hh1, hh);

        fmpz_mod_poly_get_fmpz_poly(g1, gg1);
        fmpz_poly_scalar_mul_fmpz(g1, g1, p);

        fmpz_mod_poly_get_fmpz_poly(h1, hh1);
        fmpz_poly_scalar_mul_fmpz(h1, h1, p);

        fmpz_poly_add(G, g, g1);
        fmpz_poly_add(H, h, h1);

    /* Lifting the inverses now */

        fmpz_poly_init(a1);
        fmpz_poly_init(b1);
        fmpz_poly_init(t1);
        fmpz_poly_init(t2);
        fmpz_poly_init(r);
        fmpz_poly_init(unity);

        fmpz_poly_set_coeff_si(unity, 0, -1);

        fmpz_poly_mul(t1, a, G);
        fmpz_poly_mul(t2, b, H);

        fmpz_poly_add(t1, t1, t2);
        fmpz_poly_add(t1, t1, unity);
        fmpz_poly_neg(t1, t1);

        /* Make a check that t1 is divisible by p */
        fmpz_poly_scalar_divexact_fmpz(r, t1, p);

        fmpz_mod_poly_init(rr, p1);
        fmpz_mod_poly_init(aa1, p1);
        fmpz_mod_poly_init(bb1, p1);

        fmpz_mod_poly_set_fmpz_poly(rr, r);

        /* TODO: precomputed inverses could be used here too */
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

    fmpz_poly_set(Gout, G);
    fmpz_poly_set(Hout, H);
    fmpz_poly_set(Aout, A);
    fmpz_poly_set(Bout, B);

    fmpz_poly_clear(c);
    fmpz_poly_clear(g1);
    fmpz_poly_clear(h1);
    fmpz_poly_clear(G);
    fmpz_poly_clear(H);
    fmpz_poly_clear(A);
    fmpz_poly_clear(B);

    fmpz_mod_poly_clear(cc);
    fmpz_mod_poly_clear(gg);
    fmpz_mod_poly_clear(hh);
    fmpz_mod_poly_clear(aa);
    fmpz_mod_poly_clear(bb);
    fmpz_mod_poly_clear(tt);
    fmpz_mod_poly_clear(gg1);
    fmpz_mod_poly_clear(hh1);

    fmpz_mod_poly_clear(rr);
    fmpz_mod_poly_clear(aa1);
    fmpz_mod_poly_clear(bb1);

    fmpz_poly_clear(a1);
    fmpz_poly_clear(b1);
    fmpz_poly_clear(t1);
    fmpz_poly_clear(t2);
    fmpz_poly_clear(r);
    fmpz_poly_clear(unity);
}

