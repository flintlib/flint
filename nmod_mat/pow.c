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

/*Function to implement matrix exponent
It assumes that both the matrices are square matrices and have same dimensions and mod values
It calculates the exponent by squaring each time reducing the number of multiplications to log(pow)
instead of a naive method that uses pow multiplications. */
void
nmod_mat_pow(nmod_mat_t dest, nmod_mat_t mat, ulong pow)
{
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
    else if (pow == 2)
    {
        nmod_mat_mul(dest, mat, mat);
        return;
    }

    nmod_mat_t temp1;
    nmod_mat_init_set(temp1, mat);
    if(pow == 3)
    {
        nmod_mat_mul(temp1, mat, mat);
        nmod_mat_mul(dest, temp1, mat);
        nmod_mat_clear(temp1);
        return;
    }
    nmod_mat_t temp2;
    nmod_mat_t temp3;
    nmod_mat_init_set(temp2, mat);
    nmod_mat_init_set(temp3, mat);
    nmod_mat_one(temp3);

    while(pow > 0)
    {
        if(pow%2 == 1)
        {
            nmod_mat_mul(temp2, temp3, temp1);
            nmod_mat_set(temp3, temp2);
        }
        nmod_mat_mul(temp2, temp1, temp1);
        nmod_mat_set(temp1, temp2);
        pow /= 2;
    }
    nmod_mat_set(dest, temp3);

    nmod_mat_clear(temp1);
    nmod_mat_clear(temp2);
    nmod_mat_clear(temp3);
}
