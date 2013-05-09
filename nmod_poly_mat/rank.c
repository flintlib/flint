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


len_t
nmod_poly_mat_rank(const nmod_poly_mat_t A)
{
    nmod_poly_mat_t tmp;
    nmod_poly_t den;
    len_t rank;

    if (nmod_poly_mat_is_empty(A))
        return 0;

    nmod_poly_mat_init_set(tmp, A);
    nmod_poly_init(den, nmod_poly_mat_modulus(A));
    rank = nmod_poly_mat_fflu(tmp, den, NULL, tmp, 0);
    nmod_poly_mat_clear(tmp);
    nmod_poly_clear(den);
    return rank;
}
