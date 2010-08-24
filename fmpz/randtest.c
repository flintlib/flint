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
#include <mpir.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

void
fmpz_randtest(fmpz_t f, fmpz_randstate_t state, mp_bitcnt_t bits)
{
    bits = n_randint(bits + 1);

    if (bits <= FLINT_BITS - 2)
    {
        _fmpz_demote(f);
        *f = n_randbits(bits);
        if (n_randint(2))
            *f = -*f;
    }
    else
    {
        __mpz_struct *mpz_ptr = _fmpz_promote(f);
        mpz_rrandomb(mpz_ptr, state, bits);
        if (n_randint(2))
            mpz_neg(mpz_ptr, mpz_ptr);
    }
}

void
fmpz_randtest_unsigned(fmpz_t f, fmpz_randstate_t state, mp_bitcnt_t bits)
{
    bits = n_randint(bits + 1);

    if (bits <= FLINT_BITS - 2)
    {
        _fmpz_demote(f);
        *f = n_randbits(bits);
    }
    else
    {
        __mpz_struct *mpz_ptr = _fmpz_promote(f);
        mpz_rrandomb(mpz_ptr, state, bits);
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
