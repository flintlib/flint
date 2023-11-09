/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#if defined(_WIN64) || defined(__mips64)
# include <stdint.h> /* to enable mpfr_set_sj in mpfr.h */
#endif
#include <mpfr.h>
#include "fmpz.h"

void
fmpz_get_mpfr(mpfr_t x, const fmpz_t f, mpfr_rnd_t rnd)
{
    if (!COEFF_IS_MPZ(*f))
#if defined(_WIN64) || defined(__mips64)
        mpfr_set_sj(x, *f, rnd);
#else
        mpfr_set_si(x, *f, rnd);    /* set x to small value */
#endif
    else
        mpfr_set_z(x, COEFF_TO_PTR(*f), rnd);   /* set x to large value */
}
