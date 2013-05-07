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

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "ulong_extras.h"

void
fmpz_randtest_mod(fmpz_t f, flint_rand_t state, const fmpz_t m)
{
    fmpz_t t;

    fmpz_init(t);
    fmpz_randtest_unsigned(t, state, fmpz_bits(m) + 2);
    fmpz_mod(t, t, m);

    if (n_randlimb(state) & 1UL)
    {
        fmpz_sub(t, m, t);
        fmpz_sub_ui(t, t, 1UL);
    }

    fmpz_set(f, t);
    fmpz_clear(t);
}

void
fmpz_randtest_mod_signed(fmpz_t f, flint_rand_t state, const fmpz_t m)
{
    /* Randomly generate m/2 when included in the range */
    if ((n_randlimb(state) % 32 == 1) && (fmpz_fdiv_ui(m, 2) == 0))
    {
        fmpz_fdiv_q_ui(f, m, 2UL);
    }
    else
    {
        fmpz_t t;
        fmpz_init(t);
        fmpz_tdiv_q_ui(t, m, 2UL);
        fmpz_randtest_mod(t, state, t);
        if (n_randlimb(state) & 1UL)
        {
            fmpz_neg(t, t);
        }
        fmpz_set(f, t);
        fmpz_clear(t);
    }
}
