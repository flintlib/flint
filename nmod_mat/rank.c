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

#include <stdlib.h>
#include "flint.h"
#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_mat.h"


long
nmod_mat_rank(const nmod_mat_t A)
{
    long m, n, rank;
    long * perm;
    nmod_mat_t tmp;

    m = A->r;
    n = A->c;

    if (m == 0 || n == 0)
        return 0;

    nmod_mat_init_set(tmp, A);
    perm = malloc(sizeof(long) * m);

    rank = nmod_mat_lu(perm, tmp, 0);

    free(perm);
    nmod_mat_clear(tmp);
    return rank;
}
