/*
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef NMOD_POLY_IMPL_H
#define NMOD_POLY_IMPL_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "thread_support.h"
#include "long_extras.h"
#include "mpn_extras.h"
#include "nmod_poly_factor.h"
#include "fq_nmod.h"

/* used in powmod_binexp[_preinv].c */
static __inline__ mp_limb_t
n_powmod2_mpz(mp_limb_t a, mpz_srcptr exp, mp_limb_t n, mp_limb_t ninv)
{
    if (mpz_fits_slong_p(exp))
    {
        return n_powmod2_preinv(a, flint_mpz_get_si(exp), n, ninv);
    }
    else
    {
        mpz_t t, m;
        mp_limb_t y;

        mpz_init(t);
        mpz_init(m);

        flint_mpz_set_ui(t, a);
        flint_mpz_set_ui(m, n);

        mpz_powm(t, t, exp, m);

        y = flint_mpz_get_ui(t);

        mpz_clear(t);
        mpz_clear(m);

        return y;
    }
}

#endif
