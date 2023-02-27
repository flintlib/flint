/*
    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/******************************************************************************

    Authored 2016 by Daniel S. Roche; US Government work in the public domain. 

******************************************************************************/

#include "flint.h"
#include "fmpz_mod_poly.h"
#include "fmpz.h"
#include "fmpz_vec.h"

slong _fmpz_mod_poly_minpoly_bm(fmpz* poly, 
                                const fmpz* seq, slong len, const fmpz_t p)
{
    fmpz *buf, *curpoly, *prevpoly;
    slong curlen, prevlen;
    fmpz_t disc;
    slong i, m;

    buf = _fmpz_vec_init(len + 1);
    curpoly = poly;
    _fmpz_vec_zero(curpoly, len + 1);
    prevpoly = buf;
    fmpz_init(disc);

    fmpz_one(curpoly + 0);
    curlen = 1;
    fmpz_one(prevpoly + 0);
    prevlen = 1;
    m = -1; /* m = last switching point */

    for (i = 0; i < len; ++i)
    {
        /* compute next discrepancy */
        _fmpz_vec_dot(disc, curpoly, seq + (i - curlen + 1), curlen);
        fmpz_mod(disc, disc, p);

        if (fmpz_is_zero(disc));
        else if (i - m <= curlen - prevlen)
        {
            /* quick update; no switch, curlen doesn't change. */
            slong pos = (curlen - prevlen) - (i - m);
            _fmpz_vec_scalar_addmul_fmpz(curpoly + pos, 
                prevpoly, prevlen, disc);
        }
        else
        {
            /* switching update */
            slong pos = (i - m) - (curlen - prevlen);
            _fmpz_vec_scalar_mul_fmpz(prevpoly, prevpoly, prevlen, disc);
            _fmpz_poly_add(prevpoly + pos, 
                prevpoly + pos, FLINT_MAX(0, prevlen - pos), curpoly, curlen);
            prevlen = curlen + pos;

            fmpz_sub(disc, p, disc);
            fmpz_invmod(disc, disc, p);
            _fmpz_mod_poly_scalar_mul_fmpz(curpoly, curpoly, curlen, disc, p);

            FMPZ_VEC_SWAP(curpoly, curlen, prevpoly, prevlen);
            m = i;
        }
    }

    /* make curpoly monic and copy to poly if necessary */
    fmpz_invmod(disc, curpoly + (curlen - 1), p);
    _fmpz_mod_poly_scalar_mul_fmpz(poly, curpoly, curlen, disc, p);

    _fmpz_vec_clear(buf, len + 1);
    fmpz_clear(disc);

    return curlen;
}
