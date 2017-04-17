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

slong _fmpz_sparse_index(const fmpz_sparse_t poly, const fmpz_t e)
{
    slong left = 0;
    slong right = poly->length - 1;
    while (left <= right) {
        slong mid = (left + right) / 2;
        int cmp = fmpz_cmp(e, poly->expons + mid);
        if (cmp > 0) right = mid-1;
        else if (cmp < 0) left = mid+1;
        else return mid;
    }
    return -1L - left; /* not found; return negative index. */
}
