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

#include <stdlib.h>
#include <mpir.h>
#include <mpfr.h>
#include "flint.h"
#include "math.h"
#include "arith.h"
#include "ulong_extras.h"

void
dedekind_cosine_sum_mpfr_naive(mpfr_t sum, ulong k, ulong n)
{
    mpfr_t t, u;
    long p;
    double s;
    ulong h, q;

    mpfr_set_si(sum, 0L, MPFR_RNDN);
    mpfr_init2(t, mpfr_get_prec(sum) + 5);
    mpfr_init2(u, mpfr_get_prec(sum) + 5);

    mpfr_const_pi(u, MPFR_RNDN);
    mpfr_div_ui(u, u, 6 * k, MPFR_RNDN);

    q = 6 * k;

    for (h = 0; h < k; h++)
    {
        if (n_gcd(k, h) == 1UL)
        {
            s = dedekind_sum_coprime_d(h, k);
            p = floor((s * q) + 0.5);
            p -= 12*h*n; /* XXX: 32-bit overflow */

            if (p < 0)
                p = (-p) % (12 * k);
            else
                p = p % (12 * k);

            mpfr_mul_si(t, u, p, MPFR_RNDN);
            mpfr_cos(t, t, MPFR_RNDN);
            mpfr_add(sum, sum, t, MPFR_RNDN);
        }
    }

    mpfr_clear(t);
    mpfr_clear(u);
}
