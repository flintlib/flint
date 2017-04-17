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

int fmpz_sparse_fprint_pretty(FILE * file, const fmpz_sparse_t poly, const char * x)
{
    if (poly->length == 0) fputc('0', file);
    else {
        slong i;
        for (i=0; i<poly->length; ++i) {
            if (i > 0 && fmpz_cmp_si(poly->coeffs+i, 0)>0) fputc('+', file);
            if (fmpz_is_zero(poly->expons+i)) {
                if (fmpz_fprint(file, poly->coeffs+i) <= 0) return -1;
            }
            else if (fmpz_equal_si(poly->coeffs+i, -1)) fputc('-', file);
            else if (!fmpz_is_one(poly->coeffs+i)) {
                if (fmpz_fprint(file, poly->coeffs+i) <= 0) return -1;
                fputc('*', file);
            }
            if (! fmpz_is_zero(poly->expons+i)) {
                fputs(x, file);
                if (! fmpz_is_one(poly->expons+i)) {
                    fputc('^', file);
                    if (fmpz_fprint(file, poly->expons+i) <= 0) return -1;
                }
            }
        }
    }
    return 1;
}
