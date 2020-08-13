/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "arith.h"

void mpfr_zeta_inv_euler_product(mpfr_t res, ulong s, int char_4)
{
    mpz_t z, x, y, r;
    mp_limb_t p;
    slong prec, powprec, yexp, shift;

    mpz_init(x);
    mpz_init(y);
    mpz_init(z);
    mpz_init(r);

    prec = mpfr_get_prec(res) + 32 + 2*FLINT_BIT_COUNT(s);

    flint_mpz_set_ui(z, UWORD(1));
    mpz_mul_2exp(z, z, prec);

    if (!char_4)
    {
        flint_mpz_set_ui(r, UWORD(1));
        mpz_mul_2exp(r, r, prec - s);
        mpz_sub(z, z, r);
    }

    p = UWORD(3);

    while (1)
    {
        slong i;
        powprec = prec - s*log(p)*1.4426950408889634 + 1;

        /* flint_printf("prime %wu, powprec %wd\n", p, powprec); */

        if (powprec < 5)
            break;

        flint_mpz_set_ui(x, p);
        flint_mpz_set_ui(y, UWORD(1));
        yexp = 0;

        /* Slow equivalent: flint_mpz_pow_ui(y, x, s) */
        flint_mpz_set_ui(y, p);
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

            if (s & (UWORD(1)<<i))
                flint_mpz_mul_ui(y, y, p);
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
