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

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "nmod_poly.h"


void
_fmpz_poly_CRT_ui_precomp(fmpz * res, const fmpz * poly1, long len1,
               const fmpz_t m1, mp_srcptr poly2, long len2, mp_limb_t m2,
                mp_limb_t m2inv, fmpz_t m1m2, mp_limb_t c, int sign)
{
    long i;

    for (i = 0; i < FLINT_MIN(len1, len2); i++)
    {
        _fmpz_CRT_ui_precomp(res + i, poly1 + i, m1,
                                poly2[i], m2, m2inv, m1m2, c, sign);
    }

    if (len2 > len1)
    {
        fmpz_t zero;
        fmpz_init(zero);
        for (i = len1; i < len2; i++)
        {
            _fmpz_CRT_ui_precomp(res + i, zero, m1,
                                    poly2[i], m2, m2inv, m1m2, c, sign);
        }
        fmpz_clear(zero);
    }

    for (i = len2; i < len1; i++)
    {
        _fmpz_CRT_ui_precomp(res + i, res + i, m1,
                            0, m2, m2inv, m1m2, c, sign);
    }
}

void
_fmpz_poly_CRT_ui(fmpz * res, const fmpz * poly1, long len1,
               const fmpz_t m1, mp_srcptr poly2, long len2, mp_limb_t m2,
                mp_limb_t m2inv, int sign)
{
    mp_limb_t c;
    fmpz_t m1m2;

    c = fmpz_fdiv_ui(m1, m2);
    c = n_invmod(c, m2);

    if (c == 0)
    {
        printf("Exception (_fmpz_poly_CRT_ui): m1 not invertible modulo m2.\n");
        abort();
    }

    fmpz_init(m1m2);
    fmpz_mul_ui(m1m2, m1, m2);

    _fmpz_poly_CRT_ui_precomp(res, poly1, len1, m1,
                                   poly2, len2, m2, m2inv, m1m2, c, sign);

    fmpz_clear(m1m2);
}

void
fmpz_poly_CRT_ui(fmpz_poly_t res, const fmpz_poly_t poly1,
                        const fmpz_t m1, const nmod_poly_t poly2, int sign)
{
    long len1 = poly1->length;
    long len2 = poly2->length;
    long len = FLINT_MAX(len1, len2);

    if (len == 0)
    {
        fmpz_poly_zero(res);
        return;
    }

    fmpz_poly_fit_length(res, len);

    _fmpz_poly_CRT_ui(res->coeffs, poly1->coeffs, poly1->length, m1,
        poly2->coeffs, poly2->length, poly2->mod.n, poly2->mod.ninv, sign);

    _fmpz_poly_set_length(res, len);
    _fmpz_poly_normalise(res);
}
