/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "longlong.h"
#include "mpn_extras.h"
#include "fmpzi.h"

void
fmpzi_mul(fmpzi_t res, const fmpzi_t x, const fmpzi_t y)
{
    fmpzi_struct * rp;
    fmpzi_t tmp;
    fmpz * t;
    fmpz * u;
    int xsmall, ysmall;

    const fmpz *a = fmpzi_realref(x);
    const fmpz *b = fmpzi_imagref(x);
    const fmpz *c = fmpzi_realref(y);
    const fmpz *d = fmpzi_imagref(y);

    const fmpz ca = *a;
    const fmpz cb = *b;
    const fmpz cc = *c;
    const fmpz cd = *d;

    if (x == y)
    {
        fmpzi_sqr(res, x);
        return;
    }

    xsmall = !COEFF_IS_MPZ(ca) && !COEFF_IS_MPZ(cb);
    ysmall = !COEFF_IS_MPZ(cc) && !COEFF_IS_MPZ(cd);

    if (xsmall && ysmall)
    {
        ulong thi, tlo, uhi, ulo, ahi, alo, bhi, blo;

        smul_ppmm(thi, tlo, ca, cc);
        smul_ppmm(uhi, ulo, cb, cd);
        sub_ddmmss(ahi, alo, thi, tlo, uhi, ulo);

        smul_ppmm(thi, tlo, ca, cd);
        smul_ppmm(uhi, ulo, cb, cc);
        add_ssaaaa(bhi, blo, thi, tlo, uhi, ulo);

        fmpz_set_signed_uiui(fmpzi_realref(res), ahi, alo);
        fmpz_set_signed_uiui(fmpzi_imagref(res), bhi, blo);
        return;
    }

    /* todo: detect pure real and imaginary operands */

    if (res == x || res == y)
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

    if (COEFF_IS_MPZ(ca) && COEFF_IS_MPZ(cb) &&
        COEFF_IS_MPZ(cc) && COEFF_IS_MPZ(cd))
    {
        mpz_srcptr ma = COEFF_TO_PTR(ca), mb = COEFF_TO_PTR(cb);
        mpz_srcptr mc = COEFF_TO_PTR(cc), md = COEFF_TO_PTR(cd);
        slong an = FLINT_ABS(ma->_mp_size), bn = FLINT_ABS(mb->_mp_size);
        slong cn = FLINT_ABS(mc->_mp_size), dn = FLINT_ABS(md->_mp_size);
        slong w = FLINT_MAX(an, bn) + FLINT_MAX(cn, dn) + 1;
        slong lzr, lzi;
        mpz_ptr mt, mu;

        mt = _fmpz_promote(t);
        mu = _fmpz_promote(u);

        flint_mpn_mul_complex(FLINT_MPZ_REALLOC(mt, w), &lzr,
                              FLINT_MPZ_REALLOC(mu, w), &lzi,
                              ma->_mp_d, an, ma->_mp_size < 0,
                              mb->_mp_d, bn, mb->_mp_size < 0,
                              mc->_mp_d, cn, mc->_mp_size < 0,
                              md->_mp_d, dn, md->_mp_size < 0);

        mt->_mp_size = lzr;
        mu->_mp_size = lzi;
        _fmpz_demote_val(t);
        _fmpz_demote_val(u);
        goto cleanup;
    }

    fmpz_mul(t, a, c);
    fmpz_submul(t, b, d);
    fmpz_mul(u, a, d);
    fmpz_addmul(u, b, c);

cleanup:
    if (res == x || res == y)
    {
        fmpzi_swap(res, tmp);
        fmpzi_clear(tmp);
    }
}
