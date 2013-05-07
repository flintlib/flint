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

    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

void
_fmpz_poly_pow_multinomial(fmpz * res, const fmpz * poly, long len, ulong e)
{
    long k, low, rlen;
    fmpz_t d, t;
    fmpz * P;
    
    rlen = (long) e * (len - 1L) + 1L;
    _fmpz_vec_zero(res, rlen);
    
    for (low = 0L; poly[low] == 0L; low++) ;
    if (low == 0L)
    {
        P = (fmpz *) poly;
    }
    else
    {
        P = (fmpz *) poly + low;
        len  -= low;
        res  += (long) e * low;
        rlen -= (long) e * low;
    }
    
    fmpz_init(d);
    fmpz_init(t);
    
    fmpz_pow_ui(res, P, e);
    
    for (k = 1; k < rlen; k++)
    {
        long i, u = -k;
        for (i = 1; i <= FLINT_MIN(k, len - 1); i++)
        {
            fmpz_mul(t, P + i, res + (k - i));
            u += (long) e + 1;
            if (u >= 0)
                fmpz_addmul_ui(res + k, t, (ulong) u);
            else
                fmpz_submul_ui(res + k, t, - ((ulong) u));
        }
        fmpz_add(d, d, P);
        fmpz_divexact(res + k, res + k, d);
    }
    
    fmpz_clear(d);
    fmpz_clear(t);
}

void
fmpz_poly_pow_multinomial(fmpz_poly_t res, const fmpz_poly_t poly, ulong e)
{
    const long len = poly->length;
    long rlen;

    if ((len < 2) | (e < 3UL))
    {
        if (e == 0UL)
            fmpz_poly_set_ui(res, 1);
        else if (len == 0)
            fmpz_poly_zero(res);
        else if (len == 1)
        {
            fmpz_poly_fit_length(res, 1);
            fmpz_pow_ui(res->coeffs, poly->coeffs, e);
            _fmpz_poly_set_length(res, 1);
        }
        else if (e == 1UL)
            fmpz_poly_set(res, poly);
        else  /* e == 2UL */
            fmpz_poly_sqr(res, poly);
        return;
    }
    
    rlen = (long) e * (len - 1) + 1;

    if (res != poly)
    {
        fmpz_poly_fit_length(res, rlen);
        _fmpz_poly_pow_multinomial(res->coeffs, poly->coeffs, len, e);
        _fmpz_poly_set_length(res, rlen);
    }
    else
    {
        fmpz_poly_t t;
        fmpz_poly_init2(t, rlen);
        _fmpz_poly_pow_multinomial(t->coeffs, poly->coeffs, len, e);
        _fmpz_poly_set_length(t, rlen);
        fmpz_poly_swap(res, t);
        fmpz_poly_clear(t);
    }
}
