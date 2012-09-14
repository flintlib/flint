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
    Copyright (C) 2012 Sebastian Pancratz
   
******************************************************************************/

#include "fq_poly.h"

void fq_poly_swap(fq_poly_t op1, fq_poly_t op2)
{
    if (op1 != op2)
    {
        long temp;
        fq_struct *temp_c;

        temp = op1->length;
        op1->length = op2->length;
        op2->length = temp;

        temp = op1->alloc;
        op1->alloc = op2->alloc;
        op2->alloc = temp;

        temp_c = op1->coeffs;
        op1->coeffs = op2->coeffs;
        op2->coeffs = temp_c;
    }
}
