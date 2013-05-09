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

    Copyright (C) 2010 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include "flint.h"
#include "ulong_extras.h"
#include "nmod_mat.h"
#include "nmod_vec.h"


void
nmod_mat_randrank(nmod_mat_t mat, flint_rand_t state, len_t rank)
{
    len_t i;
    mp_limb_t * diag;

    if (rank < 0 || rank > mat->r || rank > mat->c)
    {
        printf("Exception (nmod_mat_randrank). Impossible rank.\n");
        abort();
    }

    diag = _nmod_vec_init(rank);
    for (i = 0; i < rank; i++)
        diag[i] = 1 + n_randint(state, mat->mod.n - 1);

    nmod_mat_randpermdiag(mat, state, diag, rank);

    _nmod_vec_clear(diag);
}
