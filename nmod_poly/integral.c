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

/* Avoid computing every reciprocal */
#if FLINT64
#define PROD_TAKE4 65535UL
#define PROD_TAKE3 2642245UL
#define PROD_TAKE2 4294967295UL
#else
#define PROD_TAKE4 255UL
#define PROD_TAKE3 1625UL
#define PROD_TAKE2 65535UL
#endif

#define MUL3(xx,yy,zz) n_mulmod2_preinv(xx,\
    n_mulmod2_preinv(yy,zz,mod.n,mod.ninv),mod.n,mod.ninv);

void _nmod_poly_integral(mp_ptr x_int, mp_srcptr x, long len, nmod_t mod)
{
    mp_limb_t r;
    long k = len - 1;

    while (k > 0)
    {
        if (k > 3 && k < PROD_TAKE4)
        {
            r = n_invmod(k*(k-1)*(k-2)*(k-3), mod.n);
            x_int[k]   = MUL3(x[k-1], r, (k-1)*(k-2)*(k-3));
            x_int[k-1] = MUL3(x[k-2], r, k*(k-2)*(k-3));
            x_int[k-2] = MUL3(x[k-3], r, k*(k-1)*(k-3));
            x_int[k-3] = MUL3(x[k-4], r, k*(k-1)*(k-2));
            k -= 4;
        }
        else if (k > 2 && k < PROD_TAKE3)
        {
            r = n_invmod(k*(k-1)*(k-2), mod.n);
            x_int[k]   = MUL3(x[k-1], r, (k-1)*(k-2));
            x_int[k-1] = MUL3(x[k-2], r, k*(k-2));
            x_int[k-2] = MUL3(x[k-3], r, k*(k-1));
            k -= 3;
        }
        else if (k > 1 && k < PROD_TAKE2)
        {
            r = n_invmod(k*(k-1), mod.n);
            x_int[k]   = MUL3(x[k-1], r, k-1);
            x_int[k-1] = MUL3(x[k-2], r, k);
            k -= 2;
        }
        else
        {
            r = n_invmod(k, mod.n);
            x_int[k] = n_mulmod2_preinv(x[k-1], r, mod.n, mod.ninv);
            k -= 1;
        }
    }

    x_int[0] = 0UL;
}

void nmod_poly_integral(nmod_poly_t x_int, const nmod_poly_t x)
{
    nmod_poly_fit_length(x_int, x->length + 1);
    _nmod_poly_integral(x_int->coeffs, x->coeffs, x->length + 1, x->mod);
    x_int->length = x->length + 1;
	_nmod_poly_normalise(x_int);
}
