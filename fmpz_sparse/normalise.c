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

void _fmpz_sparse_normalise(fmpz_sparse_t poly)
{
    /* A modified insertion sort collects terms with the same exponent. */
    /* TODO: incorporate quicksort for larger sizes */
    slong cur = 1;
    /* nfree will hold the number of emptied terms not yet removed */
    slong nfree = 0; 
    for (; cur < poly->length; ++cur) {
        slong i = cur-nfree-1;
        /* find the position where cur belongs */
        for (; i >= 0 && fmpz_cmp(poly->expons + i, poly->expons + cur) < 0; --i);
        if (i >= 0 && fmpz_equal(poly->expons + i, poly->expons + cur)) {
            /* collision */
            fmpz_add(poly->coeffs + i, poly->coeffs + i, poly->coeffs + cur);
            fmpz_clear(poly->coeffs + cur);
            fmpz_clear(poly->expons + cur);
            ++nfree;
        }
        else if (i+1 < cur) {
            if (nfree == 0) {
                fmpz temp[2];
                temp[0] = poly->coeffs[cur];
                temp[1] = poly->expons[cur];
                _fmpz_sparse_vec_shift(poly, i+1, cur, 1);
                poly->coeffs[i+1] = temp[0];
                poly->expons[i+1] = temp[1];
            }
            else {
                _fmpz_sparse_vec_shift(poly, i+1, cur-nfree, 1);
                fmpz_swap(poly->coeffs+i+1, poly->coeffs+cur);
                fmpz_swap(poly->expons+i+1, poly->expons+cur);
            }
        }
    }
    poly->length -= nfree;
    nfree = 0;
    /* 2nd pass: remove zero terms */
    for (cur=0; cur < poly->length && ! fmpz_is_zero(poly->coeffs + cur); ++cur);
    while (cur < poly->length) {
        slong start = cur;
        for (; cur < poly->length && fmpz_is_zero(poly->coeffs + cur); ++cur);
        nfree += (cur - start);
        start = cur;
        for (; cur < poly->length && ! fmpz_is_zero(poly->coeffs + cur); ++cur);
        _fmpz_sparse_vec_shift(poly, start, cur, -nfree);
    }
    poly->length -= nfree;
}
