/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include <mpfr.h>
#include "flint.h"
#include "ulong_extras.h"

/* return ceil(log_b(2^n)) */
ulong n_clog_2exp(ulong n, ulong b)
{
    slong prec = FLINT_BITS;
    mpfr_t A, B, C;
    mpz_t Z;
    ulong r;

    if (n == 0)
        return 0;        

    if ((b & (b - 1)) == 0)
    {
        ulong log2b = 1;
        while (b > 2)
        {
            b = b >> 1;
            log2b++;
        }
        return n/log2b + (n%log2b != 0);
    }

    mpfr_init2(A, prec);
    mpfr_init2(B, prec);
    mpfr_init2(C, prec);
    mpz_init(Z);

    do {
        mpfr_set_prec(A, prec);
        mpfr_set_prec(B, prec);
        mpfr_set_prec(C, prec);

        /* compute A >= n/log2(b) */
        flint_mpz_set_ui(Z, n);
        mpfr_set_z(C, Z, MPFR_RNDA);
        flint_mpz_set_ui(Z, b);
        mpfr_set_z(B, Z, MPFR_RNDZ);
        mpfr_log2(B, B, MPFR_RNDZ);
        mpfr_div(A, C, B, MPFR_RNDA);

        mpfr_get_z(Z, A, MPFR_RNDA);
        r = flint_mpz_get_ui(Z);

        /* compute A <= n/log2(b) */
        flint_mpz_set_ui(Z, n);
        mpfr_set_z(C, Z, MPFR_RNDZ);
        flint_mpz_set_ui(Z, b);
        mpfr_set_z(B, Z, MPFR_RNDA);
        mpfr_log2(B, B, MPFR_RNDA);
        mpfr_div(A, C, B, MPFR_RNDZ);

        mpfr_get_z(Z, A, MPFR_RNDA);

        prec += FLINT_BITS;
    } while (r != flint_mpz_get_ui(Z));

    mpfr_clear(A);
    mpfr_clear(B);
    mpfr_clear(C);
    mpz_clear(Z);

    return r;
}

