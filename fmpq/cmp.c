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

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpq.h"

int
_fmpq_cmp(const fmpz_t p, const fmpz_t q, const fmpz_t r, const fmpz_t s)
{
    int s1, s2, res;
    mp_bitcnt_t bp, bq, br, bs;
    fmpz_t t, u;

    if (fmpz_equal(q, s))
        return fmpz_cmp(p, r);

    s1 = fmpz_sgn(p);
    s2 = fmpz_sgn(r);

    if (s1 != s2)
        return s1 < s2 ? -1 : 1;

    if (s1 == 0)
        return -s2;

    if (s2 == 0)
        return -s1;

    bp = fmpz_bits(p);
    bq = fmpz_bits(q);
    br = fmpz_bits(r);
    bs = fmpz_bits(s);

    if (bp + bs + 1 < br + bq)
        return -s1;

    if (bp + bs > br + bq + 1)
        return s1;

    fmpz_init(t);
    fmpz_init(u);

    fmpz_mul(t, p, s);
    fmpz_mul(u, q, r);

    res = fmpz_cmp(t, u);

    fmpz_clear(t);
    fmpz_clear(u);

    return res;
}

int
fmpq_cmp(const fmpq_t x, const fmpq_t y)
{
    return _fmpq_cmp(fmpq_numref(x), fmpq_denref(x),
                     fmpq_numref(y), fmpq_denref(y));
}

