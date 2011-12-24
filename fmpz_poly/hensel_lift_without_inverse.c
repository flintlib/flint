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
void _fmpz_poly_hensel_lift_without_inverse(fmpz *G, fmpz *H, 
    const fmpz *f, long lenF, 
    const fmpz *a, long lenA, const fmpz *b, long lenB, 
    const fmpz_t p, const fmpz_t p1)
{
}
*/
void fmpz_poly_hensel_lift_without_inverse(fmpz_poly_t Gout, fmpz_poly_t Hout, 
	const fmpz_poly_t f, const fmpz_poly_t g, const fmpz_poly_t h, 
    const fmpz_poly_t a, const fmpz_poly_t b, 
    const fmpz_t p, const fmpz_t p1)
{
    fmpz_poly_t c, g1, h1, G, H;
    fmpz_mod_poly_t cc, gg, hh, aa, bb, tt, gg1, hh1;
/*
    if (Gout == f || Gout == h || Gout == a || Gout == b ||
Hout == f || Hout == g || Hout == a || Hout == b)
{
    printf("WTF??\n"); abort();
}
*/
    fmpz_poly_init(c);
    fmpz_poly_init(g1);
    fmpz_poly_init(h1);
    fmpz_poly_init(G);
    fmpz_poly_init(H);

    fmpz_poly_mul(c, g, h);
    fmpz_poly_sub(c, f, c);
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

    fmpz_poly_set(Gout, G);
    fmpz_poly_set(Hout, H);

    fmpz_mod_poly_clear(cc);
    fmpz_mod_poly_clear(gg);
    fmpz_mod_poly_clear(hh);
    fmpz_mod_poly_clear(aa);
    fmpz_mod_poly_clear(bb);
    fmpz_mod_poly_clear(tt);
    fmpz_mod_poly_clear(gg1);
    fmpz_mod_poly_clear(hh1);

    fmpz_poly_clear(c);
    fmpz_poly_clear(g1);
    fmpz_poly_clear(h1);
    fmpz_poly_clear(G);
    fmpz_poly_clear(H);
}

