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

    Copyright (C) 2014 Ashish Kedia

******************************************************************************/

#include "flint.h"
#include "nmod_mat.h"
#include "ulong_extras.h"

void
nmod_mat_scalar_mul_add(nmod_mat_t dest, const nmod_mat_t X, const mp_limb_t b,
                        const nmod_mat_t Y)
{
    slong i, j;

    if (b == UWORD(0))
    {
        if (dest != X)
            nmod_mat_set(dest, X);
        return;
    }

    for (i = 0; i < X->r; i++)
    {
        for (j = 0; j < X->c; j++)
        {
            nmod_mat_entry(dest, i, j) =
                n_addmod(nmod_mat_entry(X, i, j),
                         n_mulmod2_preinv(nmod_mat_entry(Y, i, j), b, Y->mod.n,
                                          Y->mod.ninv), X->mod.n);
        }
    }
}
