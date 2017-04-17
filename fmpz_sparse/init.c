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

    Authored 2015 by Daniel S. Roche; US Government work in the public domain. 

******************************************************************************/

#include "fmpz_sparse.h"

static const slong FMPZ_SPARSE_INIT_LEN = 16;

void
fmpz_sparse_init(fmpz_sparse_t poly)
{
    poly->length = 0;
    poly->alloc = FMPZ_SPARSE_INIT_LEN;
    poly->coeffs = flint_calloc(poly->alloc, sizeof(fmpz));
    poly->expons = flint_calloc(poly->alloc, sizeof(fmpz));
}

void
fmpz_sparse_init2(fmpz_sparse_t poly, slong alloc)
{
  if(alloc)
  {
    poly->coeffs = (fmpz *) flint_calloc(alloc, sizeof(fmpz));
    poly->expons = (fmpz *) flint_calloc(alloc, sizeof(fmpz));
  }
  else
  {
    poly->coeffs = NULL;
    poly->expons = NULL;
  }

  poly->alloc = alloc;
  poly->length = 0;
}
