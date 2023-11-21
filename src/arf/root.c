/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <mpfr.h>
#include "arf.h"

int
arf_root(arf_ptr z, arf_srcptr x, ulong k, slong prec, arf_rnd_t rnd)
{
    mp_size_t xn, zn, val;
    mp_srcptr xptr;
    mp_ptr tmp, zptr;
    mpfr_t xf, zf;
    fmpz_t q, r;
    int inexact;

    if (k == 0)
    {
        arf_nan(z);
        return 0;
    }

    if (k == 1)
        return arf_set_round(z, x, prec, rnd);

    if (k == 2)
        return arf_sqrt(z, x, prec, rnd);

    if (arf_is_special(x))
    {
        if (arf_is_neg_inf(x))
            arf_nan(z);
        else
            arf_set(z, x);
        return 0;
    }

    if (ARF_SGNBIT(x))
    {
        arf_nan(z);
        return 0;
    }

    fmpz_init(q);
    fmpz_init(r);

    /* x = m * 2^e where e = qk + r */
    /* x^(1/k) = (m * 2^(qk+r))^(1/k)  */
    /* x^(1/k) = (m * 2^r)^(1/k) * 2^q  */
    fmpz_set_ui(r, k);
    fmpz_fdiv_qr(q, r, ARF_EXPREF(x), r);

    ARF_GET_MPN_READONLY(xptr, xn, x);
    zn = (prec + FLINT_BITS - 1) / FLINT_BITS;

    zf->_mpfr_d = tmp = flint_malloc(zn * sizeof(mp_limb_t));
    zf->_mpfr_prec = prec;
    zf->_mpfr_sign = 1;
    zf->_mpfr_exp = 0;

    xf->_mpfr_d = (mp_ptr) xptr;
    xf->_mpfr_prec = xn * FLINT_BITS;
    xf->_mpfr_sign = 1;
    xf->_mpfr_exp = fmpz_get_ui(r);

    inexact = mpfr_rootn_ui(zf, xf, k, arf_rnd_to_mpfr(rnd));
    inexact = (inexact != 0);

    val = 0;
    while (tmp[val] == 0)
        val++;

    ARF_GET_MPN_WRITE(zptr, zn - val, z);
    flint_mpn_copyi(zptr, tmp + val, zn - val);

    fmpz_add_si(ARF_EXPREF(z), q, zf->_mpfr_exp);

    flint_free(tmp);
    fmpz_clear(q);
    fmpz_clear(r);

    return inexact;
}

