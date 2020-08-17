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

void _fmpz_poly_hensel_lift(fmpz *G, fmpz *H, fmpz *A, fmpz *B, 
    const fmpz *f, slong lenF, 
    const fmpz *g, slong lenG, const fmpz *h, slong lenH, 
    const fmpz *a, slong lenA, const fmpz *b, slong lenB, 
    const fmpz_t p, const fmpz_t p1)
{
    _fmpz_poly_hensel_lift_without_inverse(G, H, f, lenF, g, lenG, h, lenH, 
        a, lenA, b, lenB, p, p1);

    _fmpz_poly_hensel_lift_only_inverse(A, B, G, lenG, H, lenH, 
        a, lenA, b, lenB, p, p1);
}

void fmpz_poly_hensel_lift(fmpz_poly_t G, fmpz_poly_t H, 
    fmpz_poly_t A, fmpz_poly_t B, 
    const fmpz_poly_t f, 
    const fmpz_poly_t g, const fmpz_poly_t h, 
    const fmpz_poly_t a, const fmpz_poly_t b, 
    const fmpz_t p, const fmpz_t p1)
{
    const slong lenG = g->length;
    const slong lenH = h->length;

    fmpz_poly_fit_length(G, lenG);
    fmpz_poly_fit_length(H, lenH);
    fmpz_poly_fit_length(A, lenH - 1);
    fmpz_poly_fit_length(B, lenG - 1);

    _fmpz_poly_hensel_lift(G->coeffs, H->coeffs, A->coeffs, B->coeffs, 
        f->coeffs, f->length, g->coeffs, g->length, h->coeffs, h->length, 
        a->coeffs, a->length, b->coeffs, b->length, p, p1);

    _fmpz_poly_set_length(G, lenG);
    _fmpz_poly_set_length(H, lenH);
    _fmpz_poly_set_length(A, lenH - 1);
    _fmpz_poly_set_length(B, lenG - 1);
    _fmpz_poly_normalise(A);
    _fmpz_poly_normalise(B);
}

