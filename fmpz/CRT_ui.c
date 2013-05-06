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

    Copyright (C) 2009 William Hart
    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

void
_fmpz_CRT_ui_precomp(fmpz_t out, const fmpz_t r1, const fmpz_t m1, ulong r2,
    ulong m2, mp_limb_t m2inv, const fmpz_t m1m2, mp_limb_t c, int sign)
{
    mp_limb_t r1mod, s;
    fmpz_t r1normal;
    fmpz_t tmp;

    fmpz_init(tmp);

    /* FIXME: assume r1 moved to [0, m1); add tests for this */
    if (fmpz_sgn(r1) < 0)
    {
        fmpz_init(r1normal);
        fmpz_add(r1normal, r1, m1);
    }
    else
    {
        *r1normal = *r1;
    }

    r1mod = fmpz_fdiv_ui(r1normal, m2);
    s = n_submod(r2, r1mod, m2);
    s = n_mulmod2_preinv(s, c, m2, m2inv);
    fmpz_mul_ui(tmp, m1, s);
    fmpz_add(tmp, tmp, r1normal);

    if (fmpz_sgn(r1) < 0)
        fmpz_clear(r1normal);

    if (sign)
    {
        fmpz_sub(out, tmp, m1m2);
        if (fmpz_cmpabs(tmp, out) <= 0)
            fmpz_set(out, tmp);
    }
    else
    {
        fmpz_set(out, tmp);
    }

    fmpz_clear(tmp);
}

void fmpz_CRT_ui(fmpz_t out, const fmpz_t r1, const fmpz_t m1,
    ulong r2, ulong m2, int sign)
{
    mp_limb_t c;
    fmpz_t m1m2;

    c = fmpz_fdiv_ui(m1, m2);
    c = n_invmod(c, m2);

    if (c == 0)
    {
        printf("Exception (fmpz_CRT_ui). m1 not invertible modulo m2.\n");
        abort();
    }

    fmpz_init(m1m2);
    fmpz_mul_ui(m1m2, m1, m2);

    _fmpz_CRT_ui_precomp(out, r1, m1, r2, m2, n_preinvert_limb(m2),
        m1m2, c, sign);

    fmpz_clear(m1m2);
}
