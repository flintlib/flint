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
#include <gmp.h>
#include <mpfr.h>
#include "flint.h"
#include "arith.h"
#include "ulong_extras.h"

void _arith_euler_number_zeta(fmpz_t res, ulong n)
{
    mpz_t r;
    mpfr_t t, z, pi;
    len_t prec, pi_prec;

    if (n % 2)
    {
        fmpz_zero(res);
        return;
    }

    if (n < SMALL_EULER_LIMIT)
    {
        fmpz_set_ui(res, euler_number_small[n / 2]);
        if (n % 4 == 2)
            fmpz_neg(res, res);
        return;
    }

    prec = arith_euler_number_size(n) + 10;
    pi_prec = prec + FLINT_BIT_COUNT(n);

    mpz_init(r);
    mpfr_init2(t, prec);
    mpfr_init2(z, prec);
    mpfr_init2(pi, pi_prec);

    mpz_fac_ui(r, n);
    mpfr_set_z(t, r, GMP_RNDN);
    mpfr_mul_2exp(t, t, n + 2, GMP_RNDN);

    /* pi^(n + 1) * L(n+1) */
    mpfr_zeta_inv_euler_product(z, n + 1, 1);
    mpfr_const_pi(pi, GMP_RNDN);
    mpfr_pow_ui(pi, pi, n + 1, GMP_RNDN);
    mpfr_mul(z, z, pi, GMP_RNDN);

    mpfr_div(t, t, z, GMP_RNDN);

    /* round */
    mpfr_round(t, t);
    mpfr_get_z(r, t, GMP_RNDN);
    fmpz_set_mpz(res, r);

    if (n % 4 == 2)
        fmpz_neg(res, res);

    mpz_clear(r);
    mpfr_clear(t);
    mpfr_clear(z);
    mpfr_clear(pi);
}
