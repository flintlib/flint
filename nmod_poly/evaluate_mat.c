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

    Copyright (C) 2010 William Hart

******************************************************************************/

#include <gmp.h>
#include <stdlib.h>
#include "flint.h"
#include "ulong_extras.h"
#include "nmod_poly.h"
#include "flint.h"
#include "nmod_mat.h"

void
nmod_mat_one_addmul(nmod_mat_t dest, nmod_mat_t mat, mp_limb_t c)
{
    slong i, j;

    if(dest == mat)
    {
        for(i = 0; i < mat->r ; i++)
        {
            nmod_mat_entry(mat, i, i) = n_addmod(nmod_mat_entry(mat, i, i),  c, mat->mod.n);
        }
        return;
    }
    for(i = 0; i < mat->r ; i++)
        for(j = 0; j < mat->c ; j++)
        {
            nmod_mat_entry(dest, i, j) = nmod_mat_entry(mat, i, j);
            if(i == j)
            {
                nmod_mat_entry(dest, i, i) = n_addmod(nmod_mat_entry(dest, i, i),  c, mat->mod.n);
            }
        }
}

void
_nmod_poly_evaluate_mat(nmod_mat_t dest, const mp_srcptr poly, slong len, const nmod_mat_t c)
{
    slong m = len-1;
    nmod_mat_t temp, point;
    int isalias = 0;

    if(dest==c)
    {
        nmod_mat_init_set(point, c);
        isalias = 1;
    }
    else
    {
        *point=*c;
    }

    nmod_mat_zero(dest);

    if (len == 0)
    {
        return;
    }

    if (len == 1 || nmod_mat_is_zero(point))
    {
        nmod_mat_one_addmul(dest, dest, poly[0]);
        return;
    }

    nmod_mat_one_addmul(dest, dest, poly[m]);
    nmod_mat_init_set(temp, point);

    for( m-- ; m >= 0 ; m--)
    {
        nmod_mat_mul(temp, dest, point);
        nmod_mat_one_addmul(dest, temp, poly[m]);
    }

    nmod_mat_clear(temp);
    if (isalias)
    {
        nmod_mat_clear(point);
    }
}

void
nmod_poly_evaluate_mat(nmod_mat_t dest, const nmod_poly_t poly, const nmod_mat_t c)
{
    _nmod_poly_evaluate_mat(dest, poly->coeffs, poly->length, c);
}
