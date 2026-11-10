/*
    Copyright (C) 2025, 2026 Marc Mezzarobba

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"

acb_poly_struct *
_acb_poly_vec_init(slong n)
{
    acb_poly_struct *vec = (acb_poly_struct *) flint_malloc(sizeof(acb_poly_struct) * n);

    for (slong i = 0; i < n; i++)
        acb_poly_init(vec + i);

    return vec;
}

void
_acb_poly_vec_clear(acb_poly_struct *vec, slong n)
{
    for (slong i = 0; i < n; i++)
        acb_poly_clear(vec + i);
    flint_free(vec);
}

void
_acb_poly_vec_zero(acb_poly_struct *dest, slong n)
{
    for (slong i = 0; i < n; i++)
        acb_poly_zero(dest + i);
}

void
_acb_poly_vec_set(acb_poly_struct *dest, const acb_poly_struct *src, slong n)
{
    for (slong i = 0; i < n; i++)
        acb_poly_set(dest + i, src + i);
}

void
_acb_poly_vec_set_block(acb_poly_struct *dest, const acb_poly_struct *src,
                        slong n, slong base, slong len)
{
    for (slong k = 0; k < n; k++)
    {
        if (base >= (src + k)->length)
        {
            acb_poly_zero(dest + k);
            continue;
        }
        slong len1 = FLINT_MIN(len, (src + k)->length - base);
        acb_poly_fit_length(dest + k, len1);
        _acb_vec_set((dest + k)->coeffs, (src + k)->coeffs + base, len1);
        _acb_poly_set_length(dest + k, len1);
        _acb_poly_normalise(dest + k);
    }
}

void
_acb_poly_vec_fit_length(acb_poly_struct *vec, slong n, slong len)
{
    for (slong i = 0; i < n; i++)
        acb_poly_fit_length(vec + i, len);
}

void
_acb_poly_vec_set_length(acb_poly_struct *vec, slong n, slong len)
{
    for (slong i = 0; i < n; i++)
        _acb_poly_set_length(vec + i, len);
}

void
_acb_poly_vec_normalise(acb_poly_struct *vec, slong n)
{
    for (slong i = 0; i < n; i++)
        _acb_poly_normalise(vec + i);
}

int
_acb_poly_vec_overlaps(const acb_poly_struct *vec1, const acb_poly_struct *vec2,
                       slong n)
{
    for (slong i = 0; i < n; i++)
        if (!acb_poly_overlaps(vec1 + i, vec2 + i))
            return 0;
    return 1;
}

slong
_acb_poly_vec_length(const acb_poly_struct * vec, slong n)
{
    slong len = 0;
    for (slong i = 0; i < n; i++)
        len = FLINT_MAX(len, vec[i].length);
    return len;
}
