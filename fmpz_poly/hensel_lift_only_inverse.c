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
#include "fmpz_mod_poly.h"

/*
    Macro for the lift B := [{(1 - aG - bH)/p} * b mod g] p + b, 
    of length at most lenG - 1.

    Assumes that {C, lenC} contains the inner part {(1 - aG - bH)/p} mod p1, 
    where lenC = max(lenA + lenG - 1, lenB + lenH - 1).  Requires temporary 
    space M, D, E.  We really only need 
        lenM = max(lenG, lenH)
        lenE = max(lenG + lenB - 2, lenH + lenA - 2)
        lenD = max(lenC, lenE)

    Writes {B, lenG - 1}.  The cofactor that is lifted is the 
    polynomial {b, lenB}, which may be aliased with B.  Although 
    it suffices to have g modulo p, there is no harm in supplying 
    {g, lenG} only reduced modulo p p1.
 */
#define liftinv(B, b, lenB, g, lenG)                                  \
do {                                                                  \
    _fmpz_vec_scalar_mod_fmpz(M, g, lenG, p1);                        \
    _fmpz_mod_poly_rem(D, C, lenC, M, lenG, one, p1);                 \
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
    _fmpz_poly_add(B, M, lenG - 1, b, lenB);                          \
} while (0)

void _fmpz_poly_hensel_lift_only_inverse(fmpz *A, fmpz *B, 
    const fmpz *G, len_t lenG, const fmpz *H, len_t lenH, 
    const fmpz *a, len_t lenA, const fmpz *b, len_t lenB, 
    const fmpz_t p, const fmpz_t p1)
{
    const fmpz one[1] = {1L};
    const len_t lenC = FLINT_MAX(lenA + lenG - 1, lenB + lenH - 1);
    const len_t lenM = FLINT_MAX(lenG, lenH);
    const len_t lenE = FLINT_MAX(lenG + lenB - 2, lenH + lenA - 2);
    const len_t lenD = FLINT_MAX(lenC, lenE);
    fmpz *C, *D, *E, *M;

    C = _fmpz_vec_init(lenC + lenD + lenD + lenM);
    D = C + lenC;
    E = D + lenD;
    M = E + lenE;

    if (lenG >= lenA)
        _fmpz_poly_mul(C, G, lenG, a, lenA);
    else
        _fmpz_poly_mul(C, a, lenA, G, lenG);
    if (lenH >= lenB)
        _fmpz_poly_mul(D, H, lenH, b, lenB);
    else
        _fmpz_poly_mul(D, b, lenB, H, lenH);
    _fmpz_vec_add(C, C, D, lenC);
    fmpz_sub_ui(C, C, 1);
    _fmpz_vec_neg(C, C, lenC);
    _fmpz_vec_scalar_divexact_fmpz(D, C, lenC, p);
    _fmpz_vec_scalar_mod_fmpz(C, D, lenC, p1);

    liftinv(B, b, lenB, G, lenG);
    liftinv(A, a, lenA, H, lenH);

    _fmpz_vec_clear(C, lenC + lenD + lenD + lenM);
}

void fmpz_poly_hensel_lift_only_inverse(fmpz_poly_t Aout, fmpz_poly_t Bout, 
    const fmpz_poly_t G, const fmpz_poly_t H, 
    const fmpz_poly_t a, const fmpz_poly_t b, 
    const fmpz_t p, const fmpz_t p1)
{
    fmpz_poly_fit_length(Aout, H->length - 1);
    fmpz_poly_fit_length(Bout, G->length - 1);

    _fmpz_poly_hensel_lift_only_inverse(Aout->coeffs, Bout->coeffs, 
        G->coeffs, G->length, H->coeffs, H->length, 
        a->coeffs, a->length, b->coeffs, b->length, p, p1);

    _fmpz_poly_set_length(Aout, H->length - 1);
    _fmpz_poly_set_length(Bout, G->length - 1);
    _fmpz_poly_normalise(Aout);
    _fmpz_poly_normalise(Bout);
}

