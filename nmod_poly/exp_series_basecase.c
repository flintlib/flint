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

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_poly.h"


void
_nmod_poly_exp_series_basecase(mp_ptr f, mp_srcptr h,
                                    long hlen, long n, nmod_t mod)
{
    long j, k;
    mp_ptr a;
    mp_limb_t s;

    f[0] = 1UL;

    a = _nmod_vec_init(FLINT_MIN(n, hlen));
    for (k = 1; k < FLINT_MIN(n, hlen); k++)
        a[k] = n_mulmod2_preinv(h[k], k, mod.n, mod.ninv);

    for (k = 1; k < n; k++)
    {
        s = n_mulmod2_preinv(a[1], f[k-1], mod.n, mod.ninv);
        for (j = 2; j < FLINT_MIN(k + 1, hlen); j++)
            NMOD_ADDMUL(s, a[j], f[k - j], mod);
        f[k] = n_mulmod2_preinv(s, n_invmod(k, mod.n), mod.n, mod.ninv);
    }

    _nmod_vec_clear(a);
}

void
nmod_poly_exp_series_basecase(nmod_poly_t f, const nmod_poly_t h, long n)
{
    long hlen;

    nmod_poly_fit_length(f, n);
    hlen = h->length;

    if (hlen > 0 && h->coeffs[0] != 0UL)
    {
        printf("Exception (nmod_poly_exp_series_basecase). Constant term != 0.\n");
        abort();
    }

    if (n <= 1 || hlen == 0)
    {
        if (n == 0)
        {
            nmod_poly_zero(f);
        }
        else
        {
            f->coeffs[0] = 1UL;
            f->length = 1;
        }
        return;
    }

    _nmod_poly_exp_series_basecase(f->coeffs, h->coeffs, hlen, n, f->mod);
    f->length = n;
    _nmod_poly_normalise(f);
}
