/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gmp.h"
#include "fmpz.h"
#ifdef LONGSLONG
# define flint_mpz_addmul_ui mpz_addmul_ui
# define flint_mpz_submul_ui mpz_submul_ui
#else
# include "gmpcompat.h"
#endif

FMPZ_INLINE
void flint_mpz_add_signed_uiui(mpz_ptr a, mpz_srcptr b, ulong c1, ulong c0)
{
    ulong d[2];
    ulong c2 = FLINT_SIGN_EXT(c1);
    mpz_t c;
    sub_ddmmss(d[1], d[0], c2^c1, c2^c0, c2, c2);
    c->_mp_d = d;
    c->_mp_alloc = 2;
    c->_mp_size = d[1] != 0 ? 2 : d[0] != 0;
    if (c2 != 0)
        c->_mp_size = -c->_mp_size;
    mpz_add(a, b, c);
}

void fmpz_addmul_si(fmpz_t f, const fmpz_t g, slong x)
{
    fmpz F, G;

    G = *g;
    if (x == 0 || G == 0)
        return;

    F = *f;
    if (F == 0)
    {
        fmpz_mul_si(f, g, x);
        return;
    }

    if (!COEFF_IS_MPZ(G))
    {
        ulong p1, p0;
        smul_ppmm(p1, p0, G, x);

        if (!COEFF_IS_MPZ(F))
        {
            ulong F1 = FLINT_SIGN_EXT(F);
            add_ssaaaa(p1, p0, p1, p0, F1, F);
            fmpz_set_signed_uiui(f, p1, p0);
        }
        else
        {
            mpz_mock_ptr pF = COEFF_TO_PTR(F);
            flint_mpz_add_signed_uiui((mpz_ptr) pF, (mpz_ptr) pF, p1, p0);
        }
    }
    else
    {
        mpz_mock_ptr pG = COEFF_TO_PTR(G);
        mpz_mock_ptr pF = _fmpz_promote_val(f);

        if (x < 0)
            flint_mpz_submul_ui((mpz_ptr) pF, (mpz_ptr) pG, -((ulong) x));
        else
            flint_mpz_addmul_ui((mpz_ptr) pF, (mpz_ptr) pG, x);

        _fmpz_demote_val(f);
    }
}
