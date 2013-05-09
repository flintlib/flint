/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2010 William Hart
    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"


static __inline__ mp_limb_t 
n_powmod2_mpz(mp_limb_t a, mpz_srcptr exp, mp_limb_t n, mp_limb_t ninv)
{
    if (mpz_fits_slong_p(exp))
    {
        return n_powmod2_preinv(a, mpz_get_si(exp), n, ninv);
    }
    else
    {
        mpz_t t, m;
        mp_limb_t y;
        mpz_init(t);
        mpz_init(m);
        mpz_set_ui(t, a);
        mpz_set_ui(m, n);
        mpz_powm(t, t, exp, m);
        y = mpz_get_ui(t);
        mpz_clear(t);
        mpz_clear(m);
        return y;
    }
}

void
_nmod_poly_powmod_mpz_binexp(mp_ptr res, mp_srcptr poly, 
                                mpz_srcptr e, mp_srcptr f,
                                len_t lenf, nmod_t mod)
{
    mp_ptr T, Q;
    len_t lenT, lenQ;
    len_t i;

    if (lenf == 2)
    {
        res[0] = n_powmod2_mpz(poly[0], e, mod.n, mod.ninv);
        return;
    }

    lenT = 2 * lenf - 3;
    lenQ = lenT - lenf + 1;

    T = _nmod_vec_init(lenT + lenQ);
    Q = T + lenT;

    _nmod_vec_set(res, poly, lenf - 1);

    for (i = mpz_sizeinbase(e, 2) - 2; i >= 0; i--)
    {
        _nmod_poly_mul(T, res, lenf - 1, res, lenf - 1, mod);
        _nmod_poly_divrem(Q, res, T, 2 * lenf - 3, f, lenf, mod);

        if (mpz_tstbit(e, i))
        {
            _nmod_poly_mul(T, res, lenf - 1, poly, lenf - 1, mod);
            _nmod_poly_divrem(Q, res, T, 2 * lenf - 3, f, lenf, mod);
        }
    }

    _nmod_vec_clear(T);
}


void
nmod_poly_powmod_mpz_binexp(nmod_poly_t res, 
                           const nmod_poly_t poly, mpz_srcptr e,
                           const nmod_poly_t f)
{
    mp_ptr p;
    len_t len = poly->length;
    len_t lenf = f->length;
    len_t trunc = lenf - 1;
    int pcopy = 0;

    if (lenf == 0)
    {
        printf("Exception (nmod_poly_powmod). Divide by zero.\n");
        abort();
    }

    if (mpz_sgn(e) < 0)
    {
        printf("Exception (nmod_poly_powmod). Negative exp not implemented.\n");
        abort();
    }

    if (len >= lenf)
    {
        nmod_poly_t t, r;
        nmod_poly_init_preinv(t, res->mod.n, res->mod.ninv);
        nmod_poly_init_preinv(r, res->mod.n, res->mod.ninv);
        nmod_poly_divrem(t, r, poly, f);
        nmod_poly_powmod_mpz_binexp(res, r, e, f);
        nmod_poly_clear(t);
        nmod_poly_clear(r);
        return;
    }

    if (mpz_fits_ulong_p(e))
    {
        ulong exp = mpz_get_ui(e);

        if (exp <= 2)
        {
            if (exp == 0UL)
            {
                nmod_poly_fit_length(res, 1);
                res->coeffs[0] = 1UL;
                res->length = 1;
            }
            else if (exp == 1UL)
            {
                nmod_poly_set(res, poly);
            }
            else
                nmod_poly_mulmod(res, poly, poly, f);
            return;
        }
    }

    if (lenf == 1 || len == 0)
    {
        nmod_poly_zero(res);
        return;
    }

    if (poly->length < trunc)
    {
        p = _nmod_vec_init(trunc);
        flint_mpn_copyi(p, poly->coeffs, poly->length);
        flint_mpn_zero(p + poly->length, trunc - poly->length);
        pcopy = 1;
    } else
        p = poly->coeffs;

    if ((res == poly && !pcopy) || (res == f))
    {
        nmod_poly_t t;
        nmod_poly_init2(t, poly->mod.n, trunc);
        _nmod_poly_powmod_mpz_binexp(t->coeffs,
            p, e, f->coeffs, lenf, poly->mod);
        nmod_poly_swap(res, t);
        nmod_poly_clear(t);
    }
    else
    {
        nmod_poly_fit_length(res, trunc);
        _nmod_poly_powmod_mpz_binexp(res->coeffs,
            p, e, f->coeffs, lenf, poly->mod);
    }

    if (pcopy)
        _nmod_vec_clear(p);

    res->length = trunc;
    _nmod_poly_normalise(res);
}
