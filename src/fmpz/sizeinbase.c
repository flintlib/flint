/*
    Copyright (C) 2010, 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include <gmp.h>
#include "long_extras.h"
#include "ulong_extras.h"
#include "fmpz.h"

size_t fmpz_sizeinbase(const fmpz_t f, int b)
{
    fmpz d = *f;
    mpz_srcptr z;
    long e;
    double s, lg;
    slong lo, hi;

    if (!COEFF_IS_MPZ(d))
        return z_sizeinbase(d, b);

    z = COEFF_TO_PTR(d);

    /* Power-of-2 bases: mpz_sizeinbase is exact. */
    if ((b & (b - 1)) == 0)
        return mpz_sizeinbase(z, b);

    if (b == 10 && mpz_size(z) == 1)
        return (size_t) n_nonzero_sizeinbase10(mpz_getlimbn(z, 0));

    /* Otherwise approximate log_b(|z|) in double precision. If its
       floor is robust to a tiny relative perturbation, return floor + 1
       in O(1). When |z| is near a power of b the check is inconclusive
       and we compare |f| against b^hi directly. */
    s = mpz_get_d_2exp(&e, z);
    lg = (log(fabs(s)) + (double) e * 0.69314718055994530942)
         * ((b == 10) ? 0.43429448190325176 : 1.0 / log((double) b));

    lo = (slong) (lg * (1.0 - 1e-12));
    hi = (slong) (lg * (1.0 + 1e-12));

    if (lo == hi)
        return (size_t) lo + 1;

    {
        fmpz_t pow;
        int cmp;

        fmpz_init(pow);
        fmpz_ui_pow_ui(pow, (ulong) b, (ulong) hi);
        cmp = fmpz_cmpabs(f, pow);
        fmpz_clear(pow);

        return (size_t) ((cmp >= 0) ? hi + 1 : hi);
    }
}
