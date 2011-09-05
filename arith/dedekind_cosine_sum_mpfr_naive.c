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
#include <math.h>
#include <mpir.h>
#include <mpfr.h>
#include "flint.h"
#include "arith.h"
#include "ulong_extras.h"

void
dedekind_cosine_sum_mpfr_naive(mpfr_t sum, mp_limb_t k, mp_limb_t n)
{
    mpfr_t t, u;
    double s, q;
    mp_limb_t h;

    mpfr_set_si(sum, 0L, MPFR_RNDN);
    mpfr_init2(t, mpfr_get_prec(sum) + 5);
    mpfr_init2(u, mpfr_get_prec(sum) + 5);

    mpfr_const_pi(u, MPFR_RNDN);
    mpfr_div_ui(u, u, 6 * k, MPFR_RNDN);

    q = fmod(12.0 * ((double) n), 12.0 * k);

    for (h = 0; h < k; h++)
    {
        if (n_gcd(k, h) == 1UL)
        {
            s = dedekind_sum_coprime_d(h, k);
            s = floor(s * (6.0 * k) + 0.5);
            s -= h*q;

            if (s < 0)
                s = fmod(-s, 12.0 * k);
            else
                s = fmod(s, 12.0 * k);

            mpfr_mul_d(t, u, s, MPFR_RNDN);
            mpfr_cos(t, t, MPFR_RNDN);
            mpfr_add(sum, sum, t, MPFR_RNDN);
        }
    }

    mpfr_clear(t);
    mpfr_clear(u);
}
