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
    Copyright (C) 2010 Fredrik Johansson

******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_mat.h"

void
nmod_mat_pow(nmod_mat_t dest,const nmod_mat_t mat,ulong pow)
{
    nmod_mat_t temp1, temp2, st, A;
    int isalias = 0;
    if (mat->r == 0)
    {
        return;
    }
    if (pow == 0)
    {
        nmod_mat_one(dest);
        return;
    }
    if (pow == 1)
    {
        nmod_mat_set(dest, mat);
        return;
    }

    if(dest == mat)
    {
        nmod_mat_init_set(A, mat);
        isalias = 1;
    }
    else
    {
        *A = *mat;
    }

    if (pow == 2)
    {
        nmod_mat_mul(dest, A, A);
        if (isalias)
        {
            nmod_mat_clear(A);
        }
        return;
    }

    nmod_mat_init(temp1, mat->r, mat->c, mat->mod.n);

    if(pow == 3)
    {
        nmod_mat_mul(temp1, A, A);
        nmod_mat_mul(dest, temp1, A);
        nmod_mat_clear(temp1);
        if(isalias)
        {
            nmod_mat_clear(A);
        }
        return;
    }

    nmod_mat_one(dest);
    nmod_mat_init_set(temp2, A);

    while(pow > 0)
    {
        if(pow%2 == 1)
        {
            nmod_mat_mul(temp1, dest, temp2);
            *st = *temp1;
            *temp1 = *dest;
            *dest = *st;
        }
        if(pow > 1)
        {
            nmod_mat_mul(temp1, temp2, temp2);
            *st = *temp1;
            *temp1 = *temp2;
            *temp2 = *st;
        }
        pow /= 2;
    }
    nmod_mat_clear(temp1);
    nmod_mat_clear(temp2);
    if(isalias)
    {
        nmod_mat_clear(A);
    }
}
