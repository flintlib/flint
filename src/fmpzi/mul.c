/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpzi.h"

void
fmpzi_mul(fmpzi_t res, const fmpzi_t x, const fmpzi_t y)
{
    slong asize, bsize, csize, dsize;
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

    if (!xsmall && !ysmall)
    {
        asize = fmpz_size(a);

        if (asize >= 13)
        {
            bsize = fmpz_size(b);
            csize = fmpz_size(c);
            dsize = fmpz_size(d);

            if (csize >= 13 && FLINT_ABS(asize - bsize) <= 2 && FLINT_ABS(csize - dsize) <= 2)
            {
                fmpz_t v;

                fmpz_init(v);
                fmpz_add(t, a, b);
                fmpz_add(v, c, d);
                fmpz_mul(u, t, v);
                fmpz_mul(t, a, c);
                fmpz_mul(v, b, d);
                fmpz_sub(u, u, t);
                fmpz_sub(u, u, v);
                fmpz_sub(t, t, v);
                fmpz_clear(v);

                goto cleanup;
            }
        }
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
