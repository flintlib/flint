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

#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

void 
_fmpz_poly_evaluate_divconquer_fmpz(fmpz_t res, const fmpz * poly, len_t len, 
                                                const fmpz_t x)
{
    len_t c, h, i, k = 1;
    fmpz *y, *T, *t = res, *u;

    h = FLINT_BIT_COUNT(len - 1);  /* 2^{h-1} < len <= 2^h */
    y = _fmpz_vec_init(2 * h + 2); /* x^{2^0}, x^{2^1}, ..., x^{2^{h-1}} */
    T = y + h;
    u = y + 2 * h + 1;

    *y = *x;
    for (i = 1; i < h; i++)
        fmpz_mul(y + i, y + (i - 1), y + (i - 1));

    for (i = 0; i < len - 1; )
    {
        fmpz_mul(u, y + 0, poly + i + 1);
        fmpz_add(t, poly + i, u);
        i += 2;
        count_trailing_zeros(c, i);
        for (k = 1; k < c; k++)
        {
            fmpz_mul(u, y + k, t);
            fmpz_add(t, T + k, u);
        }
        fmpz_swap(T + k, t);
    }
    if (len & 1L)
    {
        fmpz_set(t, poly + (len - 1));
        count_trailing_zeros(c, len + 1);
        for (k = 1; k < c; k++)
        {
            fmpz_mul(u, y + k, t);
            fmpz_add(t, T + k, u);
        }
        fmpz_swap(T + k, t);
    }
    fmpz_swap(t, T + k);

    for ( ; k < h; k++)
    {
        if ((len - 1) & (1L << k))
        {
            fmpz_mul(u, y + k, t);
            fmpz_add(t, T + k, u);
        }
    }

    *y = 0L;
    _fmpz_vec_clear(y, 2 * h + 2);
}

void
fmpz_poly_evaluate_divconquer_fmpz(fmpz_t res, const fmpz_poly_t poly, 
                                   const fmpz_t a)
{
    if (fmpz_poly_is_zero(poly))
    {
        fmpz_zero(res);
        return;
    }

    if (res == a)
    {
        fmpz_t t;

        fmpz_init(t);
        _fmpz_poly_evaluate_divconquer_fmpz(t, poly->coeffs, poly->length, a);
        fmpz_swap(res, t);
        fmpz_clear(t);
    }
    else
        _fmpz_poly_evaluate_divconquer_fmpz(res, poly->coeffs, poly->length, a);
}

