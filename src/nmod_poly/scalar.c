/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_vec.h"
#include "nmod_poly.h"

void nmod_poly_scalar_addmul_nmod(nmod_poly_t A, const nmod_poly_t B, ulong x)
{
    slong Alen = A->length;
    slong Blen = B->length;

    if (x == 0 || Blen < 1)
        return;

    nmod_poly_fit_length(A, Blen);

    if (Blen > Alen)
        _nmod_vec_zero(A->coeffs + Alen, Blen - Alen);

    _nmod_vec_scalar_addmul_nmod(A->coeffs, B->coeffs, Blen, x, A->mod);
    _nmod_poly_set_length(A, FLINT_MAX(Alen, Blen));
    _nmod_poly_normalise(A);
}

void
nmod_poly_scalar_mul_nmod(nmod_poly_t res, const nmod_poly_t poly1, ulong c)
{
    if ((poly1->length == 0) || (c == 0))
    {
        nmod_poly_zero(res);
        return;
    }

    nmod_poly_fit_length(res, poly1->length);

    _nmod_vec_scalar_mul_nmod(res->coeffs, poly1->coeffs, poly1->length,
                              c, poly1->mod);


    res->length = poly1->length;
    _nmod_poly_normalise(res);  /* there may have been cancellation */
}
