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

void _fmpz_sparse_vec_shift(fmpz_sparse_t poly, 
    slong start, slong end, slong dist)
{
    memmove(poly->coeffs+start+dist, poly->coeffs+start, (end-start)*sizeof(fmpz));
    memmove(poly->expons+start+dist, poly->expons+start, (end-start)*sizeof(fmpz));
    if (dist > 0) {
        memset(poly->coeffs+start, 0, dist*sizeof(fmpz));
        memset(poly->expons+start, 0, dist*sizeof(fmpz));
    }
    else if (dist < 0) {
        memset(poly->coeffs+end+dist, 0, (-dist)*sizeof(fmpz));
        memset(poly->expons+end+dist, 0, (-dist)*sizeof(fmpz));
    }
}
