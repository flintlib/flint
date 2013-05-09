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

    Copyright (C) 2010 William Hart
    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "nmod_poly.h"
#include "ulong_extras.h"


void
nmod_poly_deflate(nmod_poly_t result, const nmod_poly_t input, ulong deflation)
{
    len_t res_length, i;

    if (deflation == 0)
    {
        printf("Exception (nmod_poly_deflate). Division by zero.\n");
        abort();
    }

    if (input->length <= 1 || deflation == 1)
    {
        nmod_poly_set(result, input);
        return;
    }

    res_length = (input->length - 1) / deflation + 1;
    nmod_poly_fit_length(result, res_length);
    for (i = 0; i < res_length; i++)
        result->coeffs[i] = input->coeffs[i*deflation];

    result->length = res_length;
}
