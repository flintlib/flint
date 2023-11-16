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
fmpzi_sqr(fmpzi_t res, const fmpzi_t x)
{
    slong asize, bsize;
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
        asize = COEFF_TO_PTR(ca)->_mp_size;
        asize = FLINT_ABS(asize);

        if (asize >= 16)
        {
            bsize = COEFF_TO_PTR(cb)->_mp_size;
            bsize = FLINT_ABS(bsize);

            if (FLINT_ABS(asize - bsize) <= 2)
            {
                fmpz_t v;
                fmpz_init(v);

                /* a^2-b^2, (a+b)^2-a^2-b^2 */
                fmpz_add(v, a, b);
                fmpz_mul(u, v, v);
                fmpz_mul(t, a, a);
                fmpz_sub(u, u, t);
                fmpz_mul(v, b, b);
                fmpz_sub(t, t, v);
                fmpz_sub(u, u, v);

                fmpz_clear(v);

                goto cleanup;
            }
        }
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
