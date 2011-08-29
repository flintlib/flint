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

#include <math.h>
#include <mpir.h>
#include <mpfr.h>
#include "flint.h"
#include "fmpz.h"
#include "arith.h"

#define GUARD_BITS 15
#define CUTOFF 0.49
#define HEURISTIC_CUTOFF 0.01
#define HEURISTIC_COUNT 100

extern const unsigned int partitions_lookup[NUMBER_OF_SMALL_PARTITIONS];

void
number_of_partitions_heuristic(fmpz_t x, ulong n)
{
    if (n < NUMBER_OF_SMALL_PARTITIONS)
    {
        fmpz_set_ui(x, partitions_lookup[n]);
    }
    else
    {
        mpfr_t t;
        mpfr_init(t);
        number_of_partitions_mpfr(t, n, GUARD_BITS, CUTOFF, HEURISTIC_CUTOFF,
            HEURISTIC_COUNT);
        mpfr_get_z(_fmpz_promote(x), t, MPFR_RNDN);
        _fmpz_demote_val(x);
        mpfr_clear(t);
    }
}
