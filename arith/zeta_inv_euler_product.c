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
#include <math.h>
#include <gmp.h>
#include "flint.h"
#include "arith.h"
#include "ulong_extras.h"


void mpfr_zeta_inv_euler_product(mpfr_t res, ulong s, int char_4)
{
    mpz_t z, x, y, r;
    mp_limb_t p;
    len_t prec, powprec, yexp, shift;

    mpz_init(x);
    mpz_init(y);
    mpz_init(z);
    mpz_init(r);

    prec = mpfr_get_prec(res) + 32 + 2*FLINT_BIT_COUNT(s);

    mpz_set_ui(z, 1UL);
    mpz_mul_2exp(z, z, prec);

    if (!char_4)
    {
        mpz_set_ui(r, 1UL);
        mpz_mul_2exp(r, r, prec - s);
        mpz_sub(z, z, r);
    }

    p = 3UL;

    while (1)
    {
        len_t i;
        powprec = prec - s*log(p)*1.4426950408889634 + 1;

        /* printf("prime %lu, powprec %ld\n", p, powprec); */

        if (powprec < 5)
            break;

        mpz_set_ui(x, p);
        mpz_set_ui(y, 1UL);
        yexp = 0;

        /* Slow equivalent: mpz_pow_ui(y, x, s) */
        mpz_set_ui(y, p);
        for (i = FLINT_BIT_COUNT(s) - 2; i >= 0; i--)
        {
            mpz_mul(y, y, y);
            yexp += yexp;
            shift = mpz_sizeinbase(y, 2) - powprec - 4;
            if (shift >= 0)
            {
                mpz_tdiv_q_2exp(y, y, shift);
                yexp += shift;
            }

            if (s & (1UL<<i))
                mpz_mul_ui(y, y, p);
        }

        shift = yexp;
        if (shift >= 0)
            mpz_tdiv_q_2exp(r, z, shift);
        else
            mpz_mul_2exp(r, z, -shift);

        mpz_tdiv_q(r, r, y);

        if (char_4 && (p % 4 == 3))
            mpz_add(z, z, r);
        else
            mpz_sub(z, z, r);

        p = n_nextprime(p, 0);
    }

    mpfr_set_z_2exp(res, z, -prec, GMP_RNDN);

    mpz_clear(x);
    mpz_clear(y);
    mpz_clear(z);
    mpz_clear(r);
}
