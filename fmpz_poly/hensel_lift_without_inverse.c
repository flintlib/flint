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
#include "fmpz_mod_poly.h"

/*
    Macro for the lift G := [{(f - gh)/p} * b mod g] p + g.

    Assumes that {C, lenF} contains the inner part {f - gh}/p mod p1. 
    Requires temporary space M, D, E.  We really only need 
        lenM = max(lenG, lenH)
        lenE = max(lenG + lenB - 2, lenH + lenA - 2)
        lenD = max(lenE, lenF)

    Only supports aliasing between G and g.
 */
#define lift(G, g, lenG, b, lenB)                                     \
do {                                                                  \
    _fmpz_vec_scalar_mod_fmpz(M, g, lenG, p1);                        \
    _fmpz_mod_poly_rem(D, C, lenF, M, lenG, one, p1);                 \
    _fmpz_mod_poly_mul(E, D, lenG - 1, b, lenB, p1);                  \
    if (lenB > 1)                                                     \
    {                                                                 \
        _fmpz_mod_poly_rem(D, E, lenG + lenB - 2, M, lenG, one, p1);  \
        _fmpz_vec_scalar_mul_fmpz(M, D, lenG - 1, p);                 \
    }                                                                 \
    else                                                              \
    {                                                                 \
        _fmpz_vec_scalar_mul_fmpz(M, E, lenG - 1, p);                 \
    }                                                                 \
    _fmpz_vec_add(G, g, M, lenG - 1);                                 \
    fmpz_one(G + lenG - 1);                                           \
} while (0)

void _fmpz_poly_hensel_lift_without_inverse(fmpz *G, fmpz *H, 
    const fmpz *f, slong lenF, 
    const fmpz *g, slong lenG, const fmpz *h, slong lenH, 
    const fmpz *a, slong lenA, const fmpz *b, slong lenB, 
    const fmpz_t p, const fmpz_t p1)
{
    const fmpz one[1] = {1l};
    const slong lenM = FLINT_MAX(lenG, lenH);
    const slong lenE = FLINT_MAX(lenG + lenB - 2, lenH + lenA - 2);
    const slong lenD = FLINT_MAX(lenE, lenF);
    fmpz *C, *D, *E, *M;

    C = _fmpz_vec_init(lenF + lenD + lenE + lenM);
    D = C + lenF;
    E = D + lenD;
    M = E + lenE;

    if (lenG >= lenH)
        _fmpz_poly_mul(C, g,lenG, h, lenH);
    else
        _fmpz_poly_mul(C, h, lenH, g, lenG);
    _fmpz_vec_sub(C, f, C, lenF);
    _fmpz_vec_scalar_divexact_fmpz(D, C, lenF, p);
    _fmpz_vec_scalar_mod_fmpz(C, D, lenF, p1);

    lift(G, g, lenG, b, lenB);

    lift(H, h, lenH, a, lenA);

    _fmpz_vec_clear(C, lenF + lenD + lenE + lenM);
}

void fmpz_poly_hensel_lift_without_inverse(fmpz_poly_t Gout, fmpz_poly_t Hout, 
	const fmpz_poly_t f, const fmpz_poly_t g, const fmpz_poly_t h, 
    const fmpz_poly_t a, const fmpz_poly_t b, 
    const fmpz_t p, const fmpz_t p1)
{
    fmpz_poly_fit_length(Gout, g->length);
    fmpz_poly_fit_length(Hout, h->length);
    _fmpz_poly_set_length(Gout, g->length);
    _fmpz_poly_set_length(Hout, h->length);

    _fmpz_poly_hensel_lift_without_inverse(Gout->coeffs, Hout->coeffs, 
        f->coeffs, f->length, g->coeffs, g->length, h->coeffs, h->length, 
        a->coeffs, a->length, b->coeffs, b->length, p, p1);
}

