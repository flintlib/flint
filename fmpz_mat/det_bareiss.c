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
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"
#include "perm.h"

void
_fmpz_mat_det_bareiss(fmpz_t det, fmpz_mat_t tmp)
{
    len_t *perm, n = fmpz_mat_nrows(tmp);
    perm = _perm_init(n);

    fmpz_mat_fflu(tmp, det, perm, tmp, 1);

    if (_perm_parity(perm, n) == 1)
        fmpz_neg(det, det);

    _perm_clear(perm);
}

void
fmpz_mat_det_bareiss(fmpz_t det, const fmpz_mat_t A)
{
    fmpz_mat_t tmp;

    if (A->r < 1)
    {
        fmpz_one(det);
        return;
    }

    fmpz_mat_init_set(tmp, A);
    _fmpz_mat_det_bareiss(det, tmp);
    fmpz_mat_clear(tmp);
}
