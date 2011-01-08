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

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <mpir.h>

#include "flint.h"
#include "fmpz.h"
#include "ulong_extras.h"

void
fmpz_randtest(fmpz_t f, fmpz_randstate_t state, mp_bitcnt_t bits)
{
    ulong m;

    fmpz_randtest_unsigned(f, state, bits);

    m = n_randlimb();
    if (m & 1UL)
        fmpz_neg(f, f);
}

void
fmpz_randtest_unsigned(fmpz_t f, fmpz_randstate_t state, mp_bitcnt_t bits)
{
    ulong m;

    m    = n_randlimb();
    bits = n_randint(bits + 1);

    if (bits <= FLINT_BITS - 2)
    {
        _fmpz_demote(f);
        if (m & 3UL)
            *f = n_randbits(bits);
        else
        {
            m >>= 2;
            if (bits < FLINT_BITS - 2)
                *f = m & 1UL;
            else
                *f = (m & 2UL) ? COEFF_MAX : (m & 1UL);
        }
    }
    else
    {
        __mpz_struct *mpz_ptr = _fmpz_promote(f);
        mpz_rrandomb(mpz_ptr, state, bits);
        _fmpz_demote_val(f);
    }
}

void
fmpz_randtest_not_zero(fmpz_t f, fmpz_randstate_t state, mp_bitcnt_t bits)
{
    if (bits == 0)
    {
        printf("Exception: 0 passed to fmpz_randtest_not_zero\n");
        abort();
    }

    fmpz_randtest(f, state, bits);
    if (fmpz_is_zero(f))
        fmpz_set_ui(f, 1);
}
