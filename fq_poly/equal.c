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

    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2012 Andres Goens
   
******************************************************************************/

#include "fq_poly.h"

int fq_poly_equal(const fq_poly_t op1, const fq_poly_t op2)
{
    long i;

    if (op1 == op2)
        return 1;

    if (op1->length != op2->length)
        return 0;

    for (i = 0; i < op1->length; i++)
        if (!fq_equal(op1->coeffs + i, op2->coeffs + i))
            return 0;

    return 1;
}
