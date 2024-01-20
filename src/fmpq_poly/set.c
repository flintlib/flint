/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2010, 2011 Sebastian Pancratz
    Copyright (C) 2019 Tommy Hofmann

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"

void fmpq_poly_set(fmpq_poly_t poly1, const fmpq_poly_t poly2)
{
    if (poly1 != poly2)
    {
        slong i, len = poly2->length;

        fmpq_poly_fit_length(poly1, len);
        for (i = 0; i < len; i++)
            fmpz_set(poly1->coeffs + i, poly2->coeffs + i);
        _fmpq_poly_set_length(poly1, len);

        fmpz_set(poly1->den, poly2->den);
    }
}

void fmpq_poly_set_fmpq(fmpq_poly_t poly, const fmpq_t x)
{
    fmpq_poly_fit_length(poly, 1);
    fmpz_set(poly->coeffs, fmpq_numref(x));
    fmpz_set(poly->den, fmpq_denref(x));
    _fmpq_poly_set_length(poly, 1);
    _fmpq_poly_normalise(poly);
}

void fmpq_poly_set_fmpz(fmpq_poly_t poly, const fmpz_t x)
{
    fmpq_poly_fit_length(poly, 1);
    fmpz_set(poly->coeffs, x);
    fmpz_one(poly->den);
    _fmpq_poly_set_length(poly, 1);
    _fmpq_poly_normalise(poly);
}

void fmpq_poly_set_fmpz_poly(fmpq_poly_t rop, const fmpz_poly_t op)
{
    if (fmpz_poly_is_zero(op))
    {
        fmpq_poly_zero(rop);
    }
    else
    {
        fmpq_poly_fit_length(rop, fmpz_poly_length(op));
        _fmpq_poly_set_length(rop, fmpz_poly_length(op));
        _fmpz_vec_set(rop->coeffs, op->coeffs, rop->length);
        fmpz_one(rop->den);
    }
}

void
fmpq_poly_set_nmod_poly(fmpq_poly_t rop, const nmod_poly_t op)
{
    slong len = op->length;

    if (len == 0)
    {
        fmpq_poly_zero(rop);
    }
    else
    {
        slong i;
        fmpz_one(rop->den);
        fmpq_poly_fit_length(rop, len);
        for (i = 0; i < len; i++)
            fmpz_set_ui_smod(rop->coeffs + i, op->coeffs[i], op->mod.n);
        _fmpq_poly_set_length(rop, len);
    }
}

void fmpq_poly_set_si(fmpq_poly_t poly, slong x)
{
    fmpq_poly_fit_length(poly, 1);
    fmpz_set_si(poly->coeffs, x);
    fmpz_one(poly->den);
    _fmpq_poly_set_length(poly, 1);
    _fmpq_poly_normalise(poly);
}

void fmpq_poly_set_ui(fmpq_poly_t poly, ulong x)
{
    fmpq_poly_fit_length(poly, 1);
    fmpz_set_ui(poly->coeffs, x);
    fmpz_one(poly->den);
    _fmpq_poly_set_length(poly, 1);
    _fmpq_poly_normalise(poly);
}
