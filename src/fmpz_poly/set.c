/*
    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_poly.h"

void
fmpz_poly_set(fmpz_poly_t poly1, const fmpz_poly_t poly2)
{
    if (poly1 != poly2)         /* Aliasing is trivial */
    {
        slong i, len = poly2->length;

        fmpz_poly_fit_length(poly1, len);

        for (i = 0; i < len; i++)
            fmpz_set(poly1->coeffs + i, poly2->coeffs + i);

        _fmpz_poly_set_length(poly1, len);
    }
}

void
fmpz_poly_set_nmod_poly(fmpz_poly_t res, const nmod_poly_t poly)
{
    slong len = poly->length;

    if (len == 0)
    {
        fmpz_poly_zero(res);
    }
    else
    {
        slong i;
        fmpz_poly_fit_length(res, len);
        for (i = 0; i < len; i++)
            fmpz_set_ui_smod(res->coeffs + i, poly->coeffs[i], poly->mod.n);
        _fmpz_poly_set_length(res, len);
    }
}

void
fmpz_poly_set_nmod_poly_unsigned(fmpz_poly_t res, const nmod_poly_t poly)
{
    slong len = poly->length;

    if (len == 0)
    {
        fmpz_poly_zero(res);
    }
    else
    {
        slong i;
        fmpz_poly_fit_length(res, len);
        for (i = 0; i < len; i++)
            fmpz_set_ui(res->coeffs + i, poly->coeffs[i]);
        _fmpz_poly_set_length(res, len);
    }
}

void
fmpz_poly_set_fmpz(fmpz_poly_t poly, const fmpz_t c)
{
    if (fmpz_is_zero(c))
        fmpz_poly_zero(poly);
    else
    {
        fmpz_poly_fit_length(poly, 1);
        fmpz_set(poly->coeffs, c);
        _fmpz_poly_set_length(poly, 1);
    }
}

void
fmpz_poly_set_si(fmpz_poly_t poly, slong c)
{
    if (c == 0)
        fmpz_poly_zero(poly);
    else
    {
        fmpz_poly_fit_length(poly, 1);
        fmpz_set_si(poly->coeffs, c);
        _fmpz_poly_set_length(poly, 1);
    }
}

void
fmpz_poly_set_ui(fmpz_poly_t poly, ulong c)
{
    if (c == UWORD(0))
        fmpz_poly_zero(poly);
    else
    {
        fmpz_poly_fit_length(poly, 1);
        fmpz_set_ui(poly->coeffs, c);
        _fmpz_poly_set_length(poly, 1);
    }
}
