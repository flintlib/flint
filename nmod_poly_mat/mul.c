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

#include "flint.h"
#include "nmod_poly.h"
#include "nmod_poly_mat.h"

#define KS_MIN_DIM 10
#define INTERPOLATE_MIN_DIM 60
#define KS_MAX_LENGTH 128

void
nmod_poly_mat_mul(nmod_poly_mat_t C, const nmod_poly_mat_t A,
    const nmod_poly_mat_t B)
{
    len_t ar, bc, br, dim;

    ar = A->r;
    br = B->r;
    bc = B->c;

    dim = FLINT_MIN(FLINT_MIN(ar, br), bc);

    if (dim < KS_MIN_DIM)
    {
        nmod_poly_mat_mul_classical(C, A, B);
    }
    else
    {
        len_t Alen, Blen;
        mp_limb_t mod = nmod_poly_mat_modulus(A);

        Alen = nmod_poly_mat_max_length(A);
        Blen = nmod_poly_mat_max_length(B);

        if ((FLINT_BIT_COUNT(mod) > FLINT_BITS / 4)
            && (dim > INTERPOLATE_MIN_DIM + n_sqrt(FLINT_MIN(Alen, Blen)))
            && (mod >= Alen + Blen - 1) && n_is_prime(mod))
            nmod_poly_mat_mul_interpolate(C, A, B);

        else if (Alen > KS_MAX_LENGTH || Blen > KS_MAX_LENGTH)
            nmod_poly_mat_mul_classical(C, A, B);
        else
            nmod_poly_mat_mul_KS(C, A, B);
    }
}
