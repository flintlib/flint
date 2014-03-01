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
instead of a naive method that uses pow multiplications*/
void nmod_mat_pow(nmod_mat_t dest, nmod_mat_t mat,slong pow)
{
    nmod_mat_t y,temp;
    nmod_mat_init_set(y,mat);
    nmod_mat_init_set(temp,mat);
    slong i,j;

    /*Make the dest matrix as an Identity matrix*/
    nmod_mat_setiden(dest);

    while(pow>0)
    {
        if(pow%2 == 1)
        {
            nmod_mat_mul(temp,dest,y);
            nmod_mat_copy(dest,temp);
        }
        nmod_mat_mul(temp,y,y);
        nmod_mat_copy(y,temp);
        pow /= 2;
    }
}
