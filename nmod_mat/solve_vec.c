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

    Copyright (C) 2010,2011 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "flint.h"
#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_mat.h"

int
nmod_mat_solve_vec(mp_ptr x, const nmod_mat_t A, mp_srcptr b)
{
    nmod_mat_t X, B;
    int result;
    len_t i, m;

    m = A->r;

    if (m == 0)
        return 1;

    /* This is a bit of a hack. There should be a function to create
       a window into a vector */
    nmod_mat_window_init(X, A, 0, 0, m, 1);
    nmod_mat_window_init(B, A, 0, 0, m, 1);

    for (i = 0; i < m; i++) X->rows[i] = x + i;
    for (i = 0; i < m; i++) B->rows[i] = (mp_ptr) (b + i);

    result = nmod_mat_solve(X, A, B);

    nmod_mat_window_clear(X);
    nmod_mat_window_clear(B);

    return result;
}
