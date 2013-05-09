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
#include "perm.h"


mp_limb_t
_nmod_mat_det(nmod_mat_t A)
{
    mp_limb_t det;
    len_t * P;

    len_t m = A->r;
    len_t rank;
    len_t i;

    P = flint_malloc(sizeof(len_t) * m);
    rank = nmod_mat_lu(P, A, 1);

    det = 0UL;

    if (rank == m)
    {
        det = 1UL;
        for (i = 0; i < m; i++)
            det = n_mulmod2_preinv(det, nmod_mat_entry(A, i, i),
                A->mod.n, A->mod.ninv);
    }

    if (_perm_parity(P, m) == 1)
        det = nmod_neg(det, A->mod);

    flint_free(P);
    return det;
}

mp_limb_t
nmod_mat_det(const nmod_mat_t A)
{
    nmod_mat_t tmp;
    mp_limb_t det;
    len_t dim = A->r;

    if (dim != A->c)
    {
        printf("Exception (nmod_mat_det). Non-square matrix.\n");
        abort();
    }

    if (dim == 0) return 1UL;
    if (dim == 1) return nmod_mat_entry(A, 0, 0);

    nmod_mat_init_set(tmp, A);
    det = _nmod_mat_det(tmp);
    nmod_mat_clear(tmp);

    return det;
}
