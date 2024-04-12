/*
    Copyright (C) 2014 Ashish Kedia

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "nmod_poly.h"
#include "nmod_mat.h"

void
nmod_mat_one_addmul(nmod_mat_t dest, const nmod_mat_t mat, mp_limb_t c)
{
    slong i, j;

    if (dest == mat)
    {
        for (i = 0; i < mat->r; i++)
        {
            nmod_mat_entry(dest, i, i) =
                n_addmod(nmod_mat_entry(mat, i, i), c, mat->mod.n);
        }
        return;
    }
    for (i = 0; i < mat->r; i++)
        for (j = 0; j < mat->c; j++)
        {
            nmod_mat_entry(dest, i, j) = nmod_mat_entry(mat, i, j);
            if (i == j)
            {
                nmod_mat_entry(dest, i, i) =
                    n_addmod(nmod_mat_entry(dest, i, i), c, mat->mod.n);
            }
        }
}

static void
_nmod_poly_evaluate_mat_horner(nmod_mat_t dest, mp_srcptr poly, slong len, const nmod_mat_t c)
{
    slong m = len-1;
    nmod_mat_t temp;

    nmod_mat_zero(dest);

    if (len == 0)
    {
        return;
    }

    if (len == 1 || nmod_mat_is_zero(c))
    {
        nmod_mat_one_addmul(dest, dest, poly[0]);
        return;
    }

    nmod_mat_init_set(temp, c);
    nmod_mat_one_addmul(dest, dest, poly[m]);

    for( m-- ; m >= 0 ; m--)
    {
        nmod_mat_mul(temp, dest, c);
        nmod_mat_one_addmul(dest, temp, poly[m]);
    }
    nmod_mat_clear(temp);
}

void
nmod_poly_evaluate_mat_horner(nmod_mat_t dest, const nmod_poly_t poly, const nmod_mat_t c)
{
    nmod_mat_t temp;
    if (dest == c)
    {
        nmod_mat_init_set(temp, c);
        _nmod_poly_evaluate_mat_horner(dest, poly->coeffs, poly->length, temp);
        nmod_mat_clear(temp);
    }
    else
    {
        _nmod_poly_evaluate_mat_horner(dest, poly->coeffs, poly->length, c);
    }
}

void
nmod_poly_evaluate_mat_paterson_stockmeyer(nmod_mat_t dest, const nmod_poly_t poly,
                       const nmod_mat_t c)
{
    slong lim = n_sqrt(poly->length), i, j, rem, quo, curr;
    nmod_mat_t tmat, *temp;

    nmod_mat_zero(dest);

    if (poly->length == 0)
    {
        return;
    }

    if (poly->length == 1 || nmod_mat_is_zero(c))
    {
        nmod_mat_one_addmul(dest, dest, poly->coeffs[0]);
        return;
    }

    temp = flint_malloc((lim + 1) * sizeof(nmod_mat_t));
    nmod_mat_init(temp[0], c->r, c->c, c->mod.n);
    nmod_mat_one(temp[0]);
    nmod_mat_init(temp[1], c->r, c->c, c->mod.n);
    nmod_mat_set(temp[1], c);
    nmod_mat_init(tmat, c->r, c->c, c->mod.n);

    for (i = 2; i <= lim; i++)
    {
        nmod_mat_init(temp[i], c->r, c->c, c->mod.n);
        nmod_mat_mul(temp[i], temp[i - 1], c);
    }

    rem = (poly->length % lim);
    quo = poly->length / lim;
    curr = poly->length - rem - 1;

    for (i = 0; i < rem; i++)
    {
        nmod_mat_scalar_addmul_ui(dest, dest,
                                temp[i], poly->coeffs[poly->length - rem + i]);
    }

    for (i = 0; i < quo; i++)
    {
        nmod_mat_mul(tmat, dest, temp[lim]);
        nmod_mat_scalar_addmul_ui(dest, tmat, temp[lim - 1], poly->coeffs[curr--]);

        for (j = 1; j < lim; j++)
        {
            nmod_mat_scalar_addmul_ui(dest, dest, temp[lim - 1 - j],
                                                          poly->coeffs[curr--]);
        }
    }

    for (i = 0; i <= lim; i++)
    {
        nmod_mat_clear(temp[i]);
    }
    nmod_mat_clear(tmat);
    flint_free(temp);
}
