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
   
******************************************************************************/

#include <mpir.h>
#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mpoly.h"

void fmpz_mpoly_init(fmpz_mpoly_t poly, ulong vars, ulong ebits)
{
    poly->coeffs = NULL;
    poly->exps   = NULL;
    
    poly->alloc  = 0;
    poly->length = 0;
    poly->vars   = vars;
    poly->ebits  = ebits;
}

void fmpz_mpoly_init2(fmpz_mpoly_t poly, ulong alloc, ulong vars, ulong ebits)
{
    if (alloc) // allocate space for alloc
    {
        poly->coeffs = (fmpz *) calloc(alloc, sizeof(fmpz));
        poly->exps   = (fmpz *) calloc(alloc, sizeof(fmpz));
    }
    else 
    {
        poly->coeffs = NULL;
        poly->exps   = NULL;
    }

    poly->alloc  = alloc;
    poly->vars   = vars;
    poly->length = 0;
    poly->ebits  = ebits; 
}
