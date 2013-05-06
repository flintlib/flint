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

/* Small enough that a numerical computation is safe */
#define DOUBLE_CUTOFF (1UL << 21)

void
arith_dedekind_sum_coprime(fmpq_t s, const fmpz_t h, const fmpz_t k)
{
    if (fmpz_cmp_ui(k, DOUBLE_CUTOFF) < 0)
    {
        double t;

        t = arith_dedekind_sum_coprime_d(*h, *k) * (6 * (*k));

        /* Round to nearest after truncation */
        if (t > 0)
            t += 0.5;
        else
            t -= 0.5;

        fmpz_set_d(fmpq_numref(s), t);
        fmpz_set_ui(fmpq_denref(s), 6UL * (*k));
        fmpq_canonicalise(s);
    }
    else
    {
        arith_dedekind_sum_coprime_large(s, h, k);
    }
}
