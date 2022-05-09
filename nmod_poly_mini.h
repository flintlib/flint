/*
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef NMOD_POLY_MINI_H
#define NMOD_POLY_MINI_H

#ifdef NMOD_POLY_MINI_INLINES_C
#define NMOD_POLY_MINI_INLINE FLINT_DLL
#else
#define NMOD_POLY_MINI_INLINE static __inline__
#endif

#include "nmod_vec.h"

#ifdef __cplusplus
extern "C" {
#endif

/* memory management ****************************************************/

FLINT_DLL void nmod_poly_init(nmod_poly_t poly, ulong n);
FLINT_DLL void nmod_poly_init_preinv(nmod_poly_t poly, ulong n, ulong ninv);

FLINT_DLL void nmod_poly_init2(nmod_poly_t poly, ulong n, slong alloc);
FLINT_DLL void nmod_poly_init2_preinv(nmod_poly_t poly, ulong n, ulong ninv, slong alloc);

FLINT_DLL void nmod_poly_realloc(nmod_poly_t poly, slong alloc);

FLINT_DLL void nmod_poly_clear(nmod_poly_t poly);

FLINT_DLL void nmod_poly_fit_length(nmod_poly_t poly, slong alloc);

NMOD_POLY_MINI_INLINE
void nmod_poly_init_mod(nmod_poly_t poly, const nmod_t mod)
{
    poly->coeffs = NULL;
    poly->alloc = 0;
    poly->length = 0;
    poly->mod = mod;
}

NMOD_POLY_MINI_INLINE
void nmod_poly_set_mod(nmod_poly_t poly, const nmod_t mod)
{
    poly->mod = mod;
}

NMOD_POLY_MINI_INLINE
void _nmod_poly_set_length(nmod_poly_t poly, slong len)
{
    poly->length = len;
}

NMOD_POLY_MINI_INLINE
void _nmod_poly_normalise(nmod_poly_t poly)
{
    while (poly->length && (poly->coeffs[poly->length - 1] == WORD(0)))
        poly->length--;
}

/* manipulation and assignments *****************************************/

FLINT_DLL void nmod_poly_set(nmod_poly_t a, const nmod_poly_t b);

NMOD_POLY_MINI_INLINE
void nmod_poly_swap(nmod_poly_t poly1, nmod_poly_t poly2)
{
    slong t;
    ulong_ptr tp;

    t = poly1->alloc;
    poly1->alloc = poly2->alloc;
    poly2->alloc = t;

    t = poly1->length;
    poly1->length = poly2->length;
    poly2->length = t;

    tp = poly1->coeffs;
    poly1->coeffs = poly2->coeffs;
    poly2->coeffs = tp;
}

NMOD_POLY_MINI_INLINE
void nmod_poly_zero(nmod_poly_t res)
{
    res->length = 0;
}

NMOD_POLY_MINI_INLINE
void nmod_poly_one(nmod_poly_t res)
{
    nmod_poly_fit_length(res, 1);
    res->length = (res->mod.n != UWORD(1));
    res->coeffs[0] = 1;
}

NMOD_POLY_MINI_INLINE
void nmod_poly_truncate(nmod_poly_t poly, slong len)
{
    if (poly->length > len)
    {
        poly->length = len;
        _nmod_poly_normalise(poly);
    }
}

/* parameters ***********************************************************/

NMOD_POLY_MINI_INLINE
slong nmod_poly_length(const nmod_poly_t poly)
{
    return poly->length;
}

NMOD_POLY_MINI_INLINE
slong nmod_poly_degree(const nmod_poly_t poly)
{
    return poly->length - 1;
}

NMOD_POLY_MINI_INLINE
ulong nmod_poly_modulus(const nmod_poly_t poly)
{
    return poly->mod.n;
}

/* comparison ***********************************************************/

NMOD_POLY_MINI_INLINE
int nmod_poly_is_zero(const nmod_poly_t poly)
{
    return (poly->length == 0);
}

NMOD_POLY_MINI_INLINE
int nmod_poly_is_one(const nmod_poly_t poly)
{
    return (poly->mod.n == 0) || (poly->length == 1 && poly->coeffs[0] == 1);
}

NMOD_POLY_MINI_INLINE
int nmod_poly_equal(const nmod_poly_t a, const nmod_poly_t b)
{
    if (a->length != b->length)
        return 0;

    if (a != b)
        if (!_nmod_vec_equal(a->coeffs, b->coeffs, a->length))
            return 0;

    return 1;
}

/* arithmetic *****************************************************************/

FLINT_DLL void _nmod_poly_add(ulong_ptr res, ulong_srcptr poly1, slong len1, ulong_srcptr poly2, slong len2, nmod_t mod);
FLINT_DLL void nmod_poly_add(nmod_poly_t res, const nmod_poly_t poly1, const nmod_poly_t poly2);

FLINT_DLL void _nmod_poly_sub(ulong_ptr res, ulong_srcptr poly1, slong len1, ulong_srcptr poly2, slong len2, nmod_t mod);
FLINT_DLL void nmod_poly_sub(nmod_poly_t res, const nmod_poly_t poly1, const nmod_poly_t poly2);

FLINT_DLL void _nmod_poly_mul(ulong_ptr res, ulong_srcptr poly1, slong len1, ulong_srcptr poly2, slong len2, nmod_t mod);
FLINT_DLL void nmod_poly_mul(nmod_poly_t res, const nmod_poly_t poly1, const nmod_poly_t poly2);

FLINT_DLL void _nmod_poly_div(ulong_ptr Q, ulong_srcptr A, slong lenA, ulong_srcptr B, slong lenB, nmod_t mod);
FLINT_DLL void nmod_poly_div(nmod_poly_t Q, const nmod_poly_t A, const nmod_poly_t B);

FLINT_DLL void _nmod_poly_rem(ulong_ptr R, ulong_srcptr A, slong lenA, ulong_srcptr B, slong lenB, nmod_t mod);
FLINT_DLL void nmod_poly_rem(nmod_poly_t R, const nmod_poly_t A, const nmod_poly_t B);

FLINT_DLL void _nmod_poly_divrem(ulong_ptr Q, ulong_ptr R, ulong_srcptr A, slong lenA, ulong_srcptr B, slong lenB, nmod_t mod);
FLINT_DLL void nmod_poly_divrem(nmod_poly_t Q, nmod_poly_t R, const nmod_poly_t A, const nmod_poly_t B);

/* evaluation *****************************************************************/

FLINT_DLL ulong _nmod_poly_evaluate_nmod(ulong_srcptr poly, slong len, ulong c, nmod_t mod);
FLINT_DLL ulong nmod_poly_evaluate_nmod(const nmod_poly_t poly, ulong c);

/* shifting  *****************************************************************/

FLINT_DLL void _nmod_poly_shift_left(ulong_ptr res, ulong_srcptr poly, slong len, slong k);
FLINT_DLL void nmod_poly_shift_left(nmod_poly_t res, const nmod_poly_t poly, slong k);

FLINT_DLL void _nmod_poly_shift_right(ulong_ptr res, ulong_srcptr poly, slong len, slong k);
FLINT_DLL void nmod_poly_shift_right(nmod_poly_t res, const nmod_poly_t poly, slong k);

/* miscellaneous **************************************************************/

FLINT_DLL void _nmod_poly_make_monic(ulong_ptr output, ulong_srcptr input, slong len, nmod_t mod);
FLINT_DLL void nmod_poly_make_monic(nmod_poly_t output, const nmod_poly_t input);

FLINT_DLL void _nmod_poly_taylor_shift(ulong_ptr poly, ulong c, slong len, nmod_t mod);
FLINT_DLL void nmod_poly_taylor_shift(nmod_poly_t g, const nmod_poly_t f, ulong c);

#ifdef __cplusplus
}
#endif

#endif
