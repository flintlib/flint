/*
    Copyright (C) 2014 Ashish Kedia

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include <stdlib.h>
#include "flint.h"
#include "ulong_extras.h"
#include "nmod_poly.h"
#include "flint.h"
#include "nmod_mat.h"

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
