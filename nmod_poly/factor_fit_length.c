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

    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2008, 2009 William Hart

******************************************************************************/

#include <stdlib.h>
#include "nmod_poly.h"

void nmod_poly_factor_fit_length(nmod_poly_factor_t fac, len_t len)
{
    if (len > fac->alloc)
    {
        /* At least double number of allocated coeffs */
        if (len < 2 * fac->alloc)
            len = 2 * fac->alloc;
        nmod_poly_factor_realloc(fac, len);
    }
}
