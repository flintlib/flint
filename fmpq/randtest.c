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


void
_fmpq_randtest(fmpz_t num, fmpz_t den, flint_rand_t state, mp_bitcnt_t bits)
{
    ulong x = n_randlimb(state);

    fmpz_randtest(num, state, bits);
    fmpz_randtest_not_zero(den, state, bits);

    switch (x % 16)
    {
        case 0:
            fmpz_set_si(num, 1);
            break;
        case 1:
            fmpz_set_si(num, -1);
            break;
        case 2:
            fmpz_set_si(num, 2);
            break;
        case 3:
            fmpz_set_si(num, -2);
            break;
    }

    switch ((x / 16) % 16)
    {
        case 0:
            fmpz_set_si(den, 1);
            break;
        case 2:
            fmpz_set_si(den, 2);
            break;
    }

    _fmpq_canonicalise(num, den);
}

void fmpq_randtest(fmpq_t res, flint_rand_t state, mp_bitcnt_t bits)
{
    _fmpq_randtest(&res->num, &res->den, state, bits);
}

void fmpq_randtest_not_zero(fmpq_t f, flint_rand_t state, mp_bitcnt_t bits)
{
    if (bits == 0)
    {
        printf("Exception (fmpq_randtest_not_zero). bits == 0.\n");
        abort();
    }

    do {
        fmpq_randtest(f, state, bits);
    } while (fmpq_is_zero(f));
}
