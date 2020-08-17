/*
    Copyright (C) 2011 William Hart
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "ulong_extras.h"

void
_nmod_poly_div_series(mp_ptr Q, mp_srcptr A, slong Alen,
                                mp_srcptr B, slong Blen, slong n, nmod_t mod)
{
    Blen = FLINT_MIN(Blen, n);

    if (Blen < 32 || Blen < 65 * FLINT_BIT_COUNT(mod.n))
    {
        _nmod_poly_div_series_basecase(Q, A, Alen, B, Blen, n, mod);
    }
    else
    {
        mp_ptr Binv = _nmod_vec_init(n);

        _nmod_poly_inv_series(Binv, B, Blen, n, mod);
        _nmod_poly_mullow(Q, Binv, n, A, FLINT_MIN(n, Alen), n, mod);

        _nmod_vec_clear(Binv);
    }
}

void
nmod_poly_div_series(nmod_poly_t Q, const nmod_poly_t A, 
                                    const nmod_poly_t B, slong n)
{
    slong Alen, Blen;

    Blen = B->length;

    if (n == 0 || Blen == 0 || B->coeffs[0] == 0)
    {
        flint_printf("Exception (nmod_poly_div_series). Division by zero.\n");
        flint_abort();
    }

    Alen = A->length;

    if (Alen == 0)
    {
        nmod_poly_zero(Q);
        return;
    }

    if (Q != A && Q != B)
    {
        nmod_poly_fit_length(Q, n);
        _nmod_poly_div_series(Q->coeffs, A->coeffs, Alen, B->coeffs, Blen, n, Q->mod);
    }
    else
    {
        nmod_poly_t t;
        nmod_poly_init2(t, Q->mod.n, n);
        _nmod_poly_div_series(t->coeffs, A->coeffs, Alen, B->coeffs, Blen, n, Q->mod);
        nmod_poly_swap(Q, t);
        nmod_poly_clear(t);
    }

    Q->length = n;
    _nmod_poly_normalise(Q);
}

