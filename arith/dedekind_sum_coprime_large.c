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
#include "fmpq.h"
#include "ulong_extras.h"
#include "arith.h"

void
arith_dedekind_sum_coprime_large(fmpq_t s, const fmpz_t h, const fmpz_t k)
{
    fmpz_t sigma, p, pp, hh, kk, a, t;

    int sign;

    if (fmpz_cmp_ui(k, 2UL) <= 0)
    {
        fmpq_zero(s);
        return;
    }

    sign = 1;

    fmpz_init(sigma);
    fmpz_init(hh);
    fmpz_init(kk);
    fmpz_init(p);
    fmpz_init(pp);
    fmpz_init(a);
    fmpz_init(t);

    fmpz_set_ui(p, 1UL);
    fmpz_set(hh, h);
    fmpz_set(kk, k);

    while (!fmpz_is_zero(hh))
    {
        fmpz_fdiv_qr(a, t, kk, hh);

        if (sign == 1)
            fmpz_add(sigma, sigma, a);
        else
            fmpz_sub(sigma, sigma, a);

        sign = -sign;

        /* kk, hh = hh, kk mod hh */
        fmpz_swap(kk, hh);
        fmpz_swap(hh, t);

        /* p, pp = a*p + pp, p */
        fmpz_addmul(pp, a, p);
        fmpz_swap(p, pp);
    }

    if (sign < 0)
        fmpz_sub_ui(sigma, sigma, 3UL);

    /* s = (sigma + (h - p*s) / p) / 12 */
    if (sign < 0)
        fmpz_add(fmpq_numref(s), h, pp);
    else
        fmpz_sub(fmpq_numref(s), h, pp);

    fmpz_addmul(fmpq_numref(s), sigma, p);
    fmpz_mul_ui(fmpq_denref(s), p, 12UL);

    _fmpq_canonicalise(fmpq_numref(s), fmpq_denref(s));

    fmpz_clear(sigma);
    fmpz_clear(hh);
    fmpz_clear(kk);
    fmpz_clear(p);
    fmpz_clear(pp);
    fmpz_clear(a);
    fmpz_clear(t);
}
