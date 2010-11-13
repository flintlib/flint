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
#include <stdlib.h>
#include "flint.h"
#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_mat.h"


/* TODO: avoid temporary modulo reductions */
void
_nmod_mat_solve_lu_precomp(mp_limb_t * b, mp_limb_t ** const LU, long n,
    nmod_t mod)
{
    long i, j;
    mp_limb_t r;

    for (i = 0; i < n; i++)
    {
        r = 0UL;
        for (j = 0; j < i; j++)
        {
            NMOD_ADDMUL(r, LU[i][j], b[j], mod);
        }

        b[i] = nmod_sub(b[i], r, mod);
    }

    for (i = n - 1; i >= 0; i--)
    {
        r = 0UL;
        for (j = i + 1; j < n; j++)
        {
            NMOD_ADDMUL(r, LU[i][j], b[j], mod);
        }

        r = nmod_sub(b[i], r, mod);
        b[i] = n_mulmod2_preinv(r, n_invmod(LU[i][i], mod.n), mod.n, mod.ninv);
    }
}
