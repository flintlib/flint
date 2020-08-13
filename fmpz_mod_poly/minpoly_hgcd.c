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

slong _fmpz_mod_poly_minpoly_hgcd(fmpz* poly, 
                                  const fmpz* seq, slong len, const fmpz_t p)
{
    fmpz *buf, *f, *g, *A, *B;
    fmpz* M[4];
    slong lenM[4];
    slong buflen, lenA, lenB, len_poly, leng;
    int i;

    M[0] = poly;
    buflen = 7 * len + 5; /* for f, g, A, B, M[1], M[2], M[3] */
    buf = _fmpz_vec_init(buflen);
    f = buf;
    g = f + (len + 1);
    A = g + len;
    B = A + (len + 1);
    M[1] = B + len;
    M[2] = M[1] + (len + 1);
    M[3] = M[2] + (len + 1);

    /* f = x^len */
    fmpz_one(f + len);
    /* g = reversal of seq */
    for (i = 0; i < len; ++i) fmpz_set(g + i, seq + (len - i - 1));
    leng = len;
    FMPZ_VEC_NORM(g, leng);

    _fmpz_mod_poly_hgcd(M, lenM, 
            A, &lenA, B, &lenB, f, len + 1, g, leng, p);
    len_poly = lenM[0];

    /* one more step may be necessary */
    if (len_poly <= lenB)
    {
        slong quo_len = lenA - lenB + 1;
        fmpz_invmod(buf, B + (lenB - 1), p);
        _fmpz_mod_poly_divrem(M[2], M[3], A, lenA, B, lenB, buf, p); 

        if (len_poly >= quo_len) 
        {
            _fmpz_mod_poly_mul(M[3], poly, len_poly, M[2], quo_len, p);
        }
        else
        {
            _fmpz_mod_poly_mul(M[3], M[2], quo_len, poly, len_poly, p);
        }

        len_poly += quo_len - 1;
        _fmpz_mod_poly_add(poly, M[3], len_poly, M[1], lenM[1], p);
    }

    /* make poly monic */
    fmpz_invmod(buf, poly + (len_poly - 1), p);
    _fmpz_mod_poly_scalar_mul_fmpz(poly, poly, len_poly, buf, p);

    _fmpz_vec_clear(buf, buflen);

    return len_poly;
}
