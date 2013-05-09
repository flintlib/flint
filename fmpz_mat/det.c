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


void
fmpz_mat_det(fmpz_t det, const fmpz_mat_t A)
{
    len_t dim = A->r;

    if (dim != A->c)
    {
        printf("Exception (fmpz_mat_det). Non-square matrix.\n");
        abort();
    }

    if (dim < 5)
        fmpz_mat_det_cofactor(det, A);
    else if (dim < 25)
        fmpz_mat_det_bareiss(det, A);
    else if (dim < 60)
        fmpz_mat_det_modular(det, A, 1);
    else
    {
        len_t bits = fmpz_mat_max_bits(A);

        if (dim < FLINT_ABS(bits))
            fmpz_mat_det_modular(det, A, 1);
        else
            fmpz_mat_det_modular_accelerated(det, A, 1);
    }
}
