/*
    Copyright (C) 2014 Ashish Kedia

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_mat.h"

void
_nmod_mat_pow(nmod_mat_t dest, const nmod_mat_t mat, ulong pow)
{
    nmod_mat_t temp1, temp2;
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
    if (pow == 2)
    {
        nmod_mat_mul(dest, mat, mat);
        return;
    }

    nmod_mat_init(temp1, mat->r, mat->c, mat->mod.n);

    if(pow == 3)
    {
        nmod_mat_mul(temp1, mat, mat);
        nmod_mat_mul(dest, temp1, mat);
        nmod_mat_clear(temp1);
        return;
    }

    nmod_mat_one(dest);
    nmod_mat_init_set(temp2, mat);
    while(pow > 0)
    {
        if(pow%2 == 1)
        {
            nmod_mat_mul(temp1, dest, temp2);
            nmod_mat_swap_entrywise(temp1, dest);
        }
        if (pow > 1)
        {
            nmod_mat_mul(temp1, temp2, temp2);
            nmod_mat_swap_entrywise(temp1, temp2);
        }
        pow /= 2;
    }
    nmod_mat_clear(temp1);
    nmod_mat_clear(temp2);
}
void
nmod_mat_pow(nmod_mat_t dest, const nmod_mat_t mat, ulong pow)
{
    nmod_mat_t temp;
    if (mat == dest)
    {
        nmod_mat_init_set(temp, mat);
        _nmod_mat_pow(dest, temp, pow);
        nmod_mat_clear(temp);
    }
    else
    {
        _nmod_mat_pow(dest, mat, pow);
    }
}
