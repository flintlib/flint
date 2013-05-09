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

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <stdlib.h>
#include "flint.h"
#include "nmod_poly.h"
#include "nmod_poly_mat.h"

#define KS_MIN_DIM 10
#define INTERPOLATE_MIN_DIM 80
#define KS_MAX_LENGTH 128

void
nmod_poly_mat_sqr(nmod_poly_mat_t C, const nmod_poly_mat_t A)
{
    len_t dim = A->r;

    if (dim < KS_MIN_DIM)
    {
        nmod_poly_mat_sqr_classical(C, A);
    }
    else
    {
        len_t Alen;
        mp_limb_t mod = nmod_poly_mat_modulus(A);

        Alen = nmod_poly_mat_max_length(A);

        if ((FLINT_BIT_COUNT(mod) > FLINT_BITS / 4)
            && (dim > INTERPOLATE_MIN_DIM + n_sqrt(Alen))
            && (mod >= 2 * Alen - 1) && n_is_prime(mod))
            nmod_poly_mat_sqr_interpolate(C, A);

        if (Alen > KS_MAX_LENGTH)
            nmod_poly_mat_sqr_classical(C, A);
        else
            nmod_poly_mat_sqr_KS(C, A);
    }
}
