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
    Copyright (C) 2013 Mike Hansen

******************************************************************************/

#include "fq_poly.h"

void
fq_poly_inflate(fq_poly_t result, const fq_poly_t input, ulong inflation,
                const fq_ctx_t ctx)
{
    if (input->length <= 1 || inflation == 1)
    {
        fq_poly_set(result, input, ctx);
    }
    else if (inflation == 0)
    {
        fq_t v;
        fq_init(v, ctx);
        fq_one(v, ctx);
        fq_poly_evaluate_fq(v, input, v, ctx);
        fq_poly_zero(result, ctx);
        fq_poly_set_coeff(result, 0, v, ctx);
        fq_clear(v, ctx);
    }
    else
    {
        slong i, j, res_length = (input->length - 1) * inflation + 1;

        fq_poly_fit_length(result, res_length, ctx);

        for (i = input->length - 1; i > 0; i--)
        {
            fq_set(result->coeffs + (i * inflation), input->coeffs + i, ctx);
            for (j = i * inflation - 1; j > (i - 1) * inflation; j--)
                fq_zero(result->coeffs + j, ctx);
        }
        fq_set(result->coeffs, input->coeffs, ctx);
        result->length = res_length;
    }
}
