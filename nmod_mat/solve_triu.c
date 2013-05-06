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

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_mat.h"
#include "nmod_vec.h"

void
nmod_mat_solve_triu(nmod_mat_t X, const nmod_mat_t U,
                                    const nmod_mat_t B, int unit)
{
    if (B->r < NMOD_MAT_SOLVE_TRI_ROWS_CUTOFF ||
        B->c < NMOD_MAT_SOLVE_TRI_COLS_CUTOFF)
    {
        nmod_mat_solve_triu_classical(X, U, B, unit);
    }
    else
    {
        nmod_mat_solve_triu_recursive(X, U, B, unit);
    }
}
