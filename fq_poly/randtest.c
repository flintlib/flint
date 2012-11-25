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

    Copyright (C) 2009 William Hart
    Copyright (C) 2012 Andres Goens

******************************************************************************/

#include "fq_poly.h"

void fq_poly_randtest(fq_poly_t f, flint_rand_t state, 
                      long len, const fq_ctx_t ctx)
{
    long i;
    
    fq_poly_fit_length(f, len);
    for(i = 0; i < len; i++)
    {
        fq_randtest(f->coeffs + i, state, ctx);
    }
    _fq_poly_set_length(f, len);
    _fq_poly_normalise(f);
}

void fq_poly_randtest_not_zero(fq_poly_t f, flint_rand_t state, 
                               long len, const fq_ctx_t ctx)
{
    long i;

    if (len == 0)
    {
        printf("Exception (fq_poly_randtest_not_zero).  len = 0.\n");
        abort();
    }
    
    fq_poly_randtest(f, state, len, ctx);
    for(i = 0; (i < 10) && fq_poly_is_zero(f); i++)
        fq_poly_randtest(f, state, len, ctx);
    if (fq_poly_is_zero(f))
        fq_poly_one(f);
}
