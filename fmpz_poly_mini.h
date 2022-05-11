/*
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FMPZ_POLY_MINI_H
#define FMPZ_POLY_MINI_H

#ifdef FMPZ_POLY_INLINES_C
#define FMPZ_POLY_INLINE FLINT_DLL
#else
#define FMPZ_POLY_INLINE static __inline__
#endif

#include "fmpz_vec.h"

#ifdef __cplusplus
extern "C" {
#endif

/* parameters *****************************************************************/

FMPZ_POLY_INLINE
slong fmpz_poly_length(const fmpz_poly_t poly)
{
    return poly->length;
}

FMPZ_POLY_INLINE
slong fmpz_poly_degree(const fmpz_poly_t poly)
{
    return poly->length - 1;
}

/* memory management **********************************************************/

FLINT_DLL void fmpz_poly_init(fmpz_poly_t poly);
FLINT_DLL void fmpz_poly_init2(fmpz_poly_t poly, slong alloc);

FMPZ_POLY_INLINE
void fmpz_poly_clear(fmpz_poly_t poly)
{
    if (poly->coeffs)
    {
        _fmpz_vec_demote(poly->coeffs, poly->alloc);
        flint_free(poly->coeffs);
    }
}

FLINT_DLL void _fmpz_poly_normalise(fmpz_poly_t poly);

FMPZ_POLY_INLINE
void _fmpz_poly_set_length(fmpz_poly_t poly, slong newlen)
{
    if (poly->length > newlen)
    {
        _fmpz_vec_demote(poly->coeffs + newlen, poly->length - newlen);
    }
    poly->length = newlen;
}

FLINT_DLL void fmpz_poly_fit_length(fmpz_poly_t poly, slong len);

/* assignment and basic manipulation ******************************************/

FLINT_DLL void fmpz_poly_set(fmpz_poly_t poly1, const fmpz_poly_t poly2);
FLINT_DLL void fmpz_poly_set_fmpz(fmpz_poly_t poly, const fmpz_t c);
FLINT_DLL void fmpz_poly_set_ui(fmpz_poly_t poly, ulong c);
FLINT_DLL void fmpz_poly_set_si(fmpz_poly_t poly, slong c);

FLINT_DLL void fmpz_poly_swap(fmpz_poly_t poly1, fmpz_poly_t poly2);

FMPZ_POLY_INLINE
void fmpz_poly_zero(fmpz_poly_t poly)
{
   _fmpz_poly_set_length(poly, 0);
}

FMPZ_POLY_INLINE
void fmpz_poly_one(fmpz_poly_t poly)
{
    fmpz_poly_set_ui(poly, UWORD(1));
}

FLINT_DLL void fmpz_poly_set_coeff_ui(fmpz_poly_t poly, slong n, ulong x);
FLINT_DLL void fmpz_poly_set_coeff_si(fmpz_poly_t poly, slong n, slong x);

FLINT_DLL ulong fmpz_poly_get_coeff_ui(const fmpz_poly_t poly, slong n);
FLINT_DLL slong fmpz_poly_get_coeff_si(const fmpz_poly_t poly, slong n);

#define fmpz_poly_lead(poly) \
    ((poly)->length ? (poly)->coeffs + (poly)->length - 1 : NULL)

/* comparison *****************************************************************/

FLINT_DLL int fmpz_poly_equal(const fmpz_poly_t poly1, const fmpz_poly_t poly2);

#define fmpz_poly_is_zero(poly) \
    ((poly)->length == 0)

FMPZ_POLY_INLINE
int fmpz_poly_is_one(const fmpz_poly_t op)
{
    return (op->length) == 1 && (*(op->coeffs) == WORD(1));
}

/* operations *****************************************************************/

FLINT_DLL void fmpz_poly_neg(fmpz_poly_t res, const fmpz_poly_t poly);

FLINT_DLL void _fmpz_poly_add(fmpz * res, const fmpz * poly1, slong len1, const fmpz * poly2, slong len2);
FLINT_DLL void fmpz_poly_add(fmpz_poly_t res, const fmpz_poly_t poly1, const fmpz_poly_t poly2);

FLINT_DLL void _fmpz_poly_sub(fmpz * res, const fmpz * poly1, slong len1, const fmpz * poly2, slong len2);
FLINT_DLL void fmpz_poly_sub(fmpz_poly_t res, const fmpz_poly_t poly1, const fmpz_poly_t poly2);

FLINT_DLL void _fmpz_poly_mul(fmpz * res, const fmpz * poly1, slong len1, const fmpz * poly2, slong len2);
FLINT_DLL void fmpz_poly_mul(fmpz_poly_t res, const fmpz_poly_t poly1, const fmpz_poly_t poly2);

FLINT_DLL int _fmpz_poly_div(fmpz * Q, const fmpz * A, slong lenA, const fmpz * B, slong lenB, int exact);
FLINT_DLL void fmpz_poly_div(fmpz_poly_t Q, const fmpz_poly_t A, const fmpz_poly_t B);

/* miscellaneous **************************************************************/

FMPZ_POLY_INLINE
void fmpz_poly_content(fmpz_t res, const fmpz_poly_t poly)
{
    _fmpz_vec_content(res, poly->coeffs, poly->length);
}

#ifdef __cplusplus
}
#endif

#endif
