/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "longlong.h"
#include "mpn_extras.h"
#include "fmpzi.h"

void
fmpzi_sqr(fmpzi_t res, const fmpzi_t x)
{
    fmpzi_struct * rp;
    fmpzi_t tmp;
    fmpz * t;
    fmpz * u;

    const fmpz *a = fmpzi_realref(x);
    const fmpz *b = fmpzi_imagref(x);

    const fmpz ca = *a;
    const fmpz cb = *b;

    if (!COEFF_IS_MPZ(ca) && !COEFF_IS_MPZ(cb))
    {
        ulong thi, tlo, uhi, ulo, ahi, alo, bhi, blo;

        smul_ppmm(thi, tlo, ca, ca);
        smul_ppmm(uhi, ulo, cb, cb);
        sub_ddmmss(ahi, alo, thi, tlo, uhi, ulo);
        smul_ppmm(bhi, blo, ca + ca, cb);

        fmpz_set_signed_uiui(fmpzi_realref(res), ahi, alo);
        fmpz_set_signed_uiui(fmpzi_imagref(res), bhi, blo);

        return;
    }

    if (cb == 0)
    {
        fmpz_mul(fmpzi_realref(res), a, a);
        fmpz_zero(fmpzi_imagref(res));
        return;
    }

    if (ca == 0)
    {
        fmpz_mul(fmpzi_realref(res), b, b);
        fmpz_neg(fmpzi_realref(res), fmpzi_realref(res));
        fmpz_zero(fmpzi_imagref(res));
        return;
    }

    if (res == x)
    {
        fmpzi_init(tmp);
        rp = tmp;
    }
    else
    {
        rp = res;
    }

    t = fmpzi_realref(rp);
    u = fmpzi_imagref(rp);

    if (COEFF_IS_MPZ(ca) && COEFF_IS_MPZ(cb))
    {
        mpz_srcptr ma = COEFF_TO_PTR(ca), mb = COEFF_TO_PTR(cb);
        slong an = FLINT_ABS(ma->_mp_size), bn = FLINT_ABS(mb->_mp_size);
        slong w = 2 * FLINT_MAX(an, bn) + 1;
        slong lzr, lzi;
        mpz_ptr mt, mu;

        mt = _fmpz_promote(t);
        mu = _fmpz_promote(u);

        flint_mpn_sqr_complex(FLINT_MPZ_REALLOC(mt, w), &lzr,
                              FLINT_MPZ_REALLOC(mu, w), &lzi,
                              ma->_mp_d, an, ma->_mp_size < 0,
                              mb->_mp_d, bn, mb->_mp_size < 0);

        mt->_mp_size = lzr;
        mu->_mp_size = lzi;
        _fmpz_demote_val(t);
        _fmpz_demote_val(u);
        goto cleanup;
    }

    fmpz_mul(t, a, a);
    fmpz_mul(u, b, b);
    fmpz_sub(t, t, u);
    fmpz_mul(u, a, b);
    fmpz_mul_2exp(u, u, 1);

cleanup:
    if (res == x)
    {
        fmpzi_swap(res, tmp);
        fmpzi_clear(tmp);
    }
}
