/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2009 Andy Novocin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

extern double __gmpn_get_d(mp_limb_t *, size_t, size_t, long);

double
fmpz_get_d_2exp(slong *exp, const fmpz_t f)
{
    fmpz d = *f;

    if (!COEFF_IS_MPZ(d))
    {
        ulong d_abs;
        if (d == WORD(0))
        {
            (*exp) = WORD(0);
            return 0.0;
        }
        d_abs = FLINT_ABS(d);
        (*exp) = FLINT_BIT_COUNT(d_abs);
        if (d < WORD(0))
            return __gmpn_get_d((mp_limb_t *) &d_abs, WORD(1), WORD(-1), -*exp);
        else
            return __gmpn_get_d((mp_limb_t *) &d, WORD(1), WORD(1), -*exp);
    }
    else
    {
#if defined(__MPIR_VERSION)
       return mpz_get_d_2exp(exp, COEFF_TO_PTR(d));
#else
       long exp2;
       double m = mpz_get_d_2exp(&exp2, COEFF_TO_PTR(d));
       *exp = exp2;
       return m;
#endif
    }
}
