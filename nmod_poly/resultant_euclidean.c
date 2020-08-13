/*
    Copyright (C) 2007, 2008 William Hart
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "mpn_extras.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

mp_limb_t 
_nmod_poly_resultant_euclidean(mp_srcptr poly1, slong len1, 
                               mp_srcptr poly2, slong len2, nmod_t mod)
{
    if (poly1 == poly2)
    {
        return 0;
    }
    else if (len2 == 1)
    {
        if (len1 == 1)
        {
            return 1;
        }
        else if (len1 == 2)
        {
            return poly2[0];
        }
        else
        {
            return n_powmod2_ui_preinv(poly2[0], len1 - 1, mod.n, mod.ninv);
        }
    }
    else  /* len1 >= len2 >= 2 */
    {
        mp_limb_t res = 1;

        mp_ptr u, v, r, t, w;
        slong l0, l1, l2;
        mp_limb_t lc;

        w = _nmod_vec_init(3 * len1);
        u = w;
        v = w + len1;
        r = v + len1;

        _nmod_vec_set(u, poly1, len1);
        _nmod_vec_set(v, poly2, len2);
        l1 = len1;
        l2 = len2;

        do
        {
            l0 = l1;
            l1 = l2;
            lc = v[l1 - 1];

            _nmod_poly_rem(r, u, l0, v, l1, mod);
            l2 = l1 - 1;
            MPN_NORM(r, l2);
            {
                t = u;
                u = v;
                v = r;
                r = t;
            }

            if (l2 >= 1) 
            {
                lc  = n_powmod2_preinv(lc, l0 - l2, mod.n, mod.ninv);
                res = n_mulmod2_preinv(res, lc, mod.n, mod.ninv);

                if (((l0 | l1) & 1) == 0)
                {
                    res = nmod_neg(res, mod);
                }  
            }
            else 
            {
                if (l1 == 1)
                {
                    lc  = n_powmod2_preinv(lc, l0 - 1, mod.n, mod.ninv);
                    res = n_mulmod2_preinv(res, lc, mod.n, mod.ninv);
                }
                else
                {
                    res = 0;
                }
            }
        }
        while (l2 > 0);

        _nmod_vec_clear(w);

        return res;
    }
}

mp_limb_t 
nmod_poly_resultant_euclidean(const nmod_poly_t f, const nmod_poly_t g)
{
    const slong len1 = f->length;
    const slong len2 = g->length;
    mp_limb_t r;

    if (len1 == 0 || len2 == 0)
    {
        r = 0;
    }
    else
    {
        if (len1 >= len2)
        {
            r = _nmod_poly_resultant_euclidean(f->coeffs, len1, 
                                               g->coeffs, len2, f->mod);
        }
        else
        {
            r = _nmod_poly_resultant_euclidean(g->coeffs, len2, 
                                               f->coeffs, len1, f->mod);

            if (((len1 | len2) & WORD(1)) == WORD(0))
                r = nmod_neg(r, f->mod);
        }
    }

    return r;
}

