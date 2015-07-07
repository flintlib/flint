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

    Copyright (C) 2015 Vladimir Glazachev

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

void
fmpz_mul_mont(fmpz_t f, const fmpz_t g, const fmpz_t h,
        const fmpz_t n)
{
    ulong u, g_bit, h_bit;
    slong i, len;

    len = fmpz_bits(n);
    fmpz_set_ui(f, 0);
    h_bit = fmpz_tstbit(h, 0);

    for (i = 0; i <= len - 1; i++)
    {
        g_bit = fmpz_tstbit(g, i);

        u = fmpz_tstbit(f, 0) ^ (g_bit & h_bit);

        fmpz_addmul_ui(f, h, g_bit);
        fmpz_addmul_ui(f, n, u);
        fmpz_fdiv_q_2exp(f, f, 1);
    }

    if (fmpz_cmp(f, n) >= 0)
        fmpz_sub(f, f, n);
}

