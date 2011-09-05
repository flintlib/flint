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
#include <string.h>
#include <mpir.h>
#include <mpfr.h>
#include "flint.h"
#include "arith.h"
#include "ulong_extras.h"

void
dedekind_cosine_sum_mpfr_fast(mpfr_t sum, mp_limb_t k, mp_limb_t n)
{
    trig_prod_t prod;
    mpfr_prec_t wp;
    mpfr_t t, pi;
    mp_limb_t v;
    int i;

    if (k <= 2)
    {
        if (k == 0)
            mpfr_set_ui(sum, 0UL, MPFR_RNDN);
        else if (k == 1)
            mpfr_set_ui(sum, 1UL, MPFR_RNDN);
        else if (n % 2 == 1)
            mpfr_set_si(sum, -1L, MPFR_RNDN);
        else
            mpfr_set_si(sum, 1L, MPFR_RNDN);
        return;
    }

    trig_prod_init(prod);
    dedekind_cosine_sum_factored(prod, k, n);

    if (prod->prefactor == 0)
    {
        mpfr_set_ui(sum, 0UL, MPFR_RNDN);
        return;
    }

    wp = mpfr_get_prec(sum) + 5;
    mpfr_init2(t, wp);
    mpfr_init2(pi, wp);

    mpfr_const_pi(pi, MPFR_RNDN);

    mpfr_set_si(sum, prod->prefactor, MPFR_RNDN);
    v = n_gcd(FLINT_MAX(prod->sqrt_p, prod->sqrt_q),
              FLINT_MIN(prod->sqrt_p, prod->sqrt_q));
    prod->sqrt_p /= v;
    prod->sqrt_q /= v;

    if (prod->sqrt_p != 1)
    {
        mpfr_sqrt_ui(t, prod->sqrt_p, MPFR_RNDN);
        mpfr_mul(sum, sum, t, MPFR_RNDN);
    }

    if (prod->sqrt_q != 1)
    {
        mpfr_sqrt_ui(t, prod->sqrt_q, MPFR_RNDN);
        mpfr_div(sum, sum, t, MPFR_RNDN);
    }

    for (i = 0; i < prod->n; i++)
    {
        mpfr_mul_si(t, pi, prod->cos_p[i], MPFR_RNDN);
        mpfr_div_ui(t, t, prod->cos_q[i], MPFR_RNDN);
        mpfr_cos(t, t, MPFR_RNDN);
        mpfr_mul(sum, sum, t, MPFR_RNDN);
    }

    mpfr_clear(pi);
    mpfr_clear(t);
}
