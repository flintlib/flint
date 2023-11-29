/*
    Copyright (C) 2020 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

/*
    compute f^0, f^1, ..., f^(n-1) mod g, where g has length glen and f is
    reduced mod g and has length flen (possibly zero spaced)
    assumes res is an array of n arrays each with space for at least glen - 1
    coefficients and that flen > 0
    {ginv, ginvlen} must be set to the power series inverse of the reverse of g
*/
void
_nmod_poly_powers_mod_preinv_naive(mp_ptr * res, mp_srcptr f, slong flen, slong n,
         mp_srcptr g, slong glen, mp_srcptr ginv, slong ginvlen, const nmod_t mod)
{
    slong i;

    if (n == 0)
        return;

    /* f^0 = 1 */
    if (glen > 1)
        res[0][0] = 1;

    if (glen > 2)
       flint_mpn_zero(res[0] + 1, glen - 2);

    if (n == 1)
       return;

    /* f^1 = f */
    _nmod_vec_set(res[1], f, flen);
    flint_mpn_zero(res[1] + flen, glen - flen - 1);

    if (n == 2)
       return;

    /* f^i = f^(i - 1)*f */
    if (glen == 2) /* special case, constant polys */
    {
        for (i = 2; i < n; i++)
            res[i][0] = n_mulmod2_preinv(res[i - 1][0], res[1][0],
                                                              mod.n, mod.ninv);
    } else
    {
        for (i = 2; i < n; i++)
            _nmod_poly_mulmod_preinv(res[i], res[i - 1], glen - 1, res[1],
                                        glen - 1, g, glen, ginv, ginvlen, mod);
    }
}

void
nmod_poly_powers_mod_naive(nmod_poly_struct * res, const nmod_poly_t f,
               slong n, const nmod_poly_t g)
{
    slong i;

    nmod_poly_t ginv;
    mp_ptr * res_arr;

    if (nmod_poly_length(g) == 0)
    {
        flint_throw(FLINT_ERROR, "Exception (nmod_poly_powers_mod_naive). Divide by zero.\n");
    }

    if (nmod_poly_length(f) == 0 || nmod_poly_length(g) == 1)
    {
        if (n > 0)
           nmod_poly_one(res + 0);

        for (i = 1; i < n; i++)
           nmod_poly_zero(res + i);

        return;
    }

    if (nmod_poly_length(f) >= nmod_poly_length(g))
    {
        nmod_poly_t q, r;

        nmod_poly_init_mod(q, f->mod);
        nmod_poly_init_mod(r, f->mod);

        nmod_poly_divrem(q, r, f, g);
        nmod_poly_powers_mod_naive(res, r, n, g);

        nmod_poly_clear(q);
        nmod_poly_clear(r);

        return;
    }

    res_arr = (mp_ptr *) flint_malloc(n*sizeof(mp_ptr));
    nmod_poly_init_mod(ginv, g->mod);

    for (i = 0; i < n; i++)
    {
       nmod_poly_fit_length(res + i, nmod_poly_length(g) - 1);
       res_arr[i] = res[i].coeffs;
       _nmod_poly_set_length(res + i, nmod_poly_length(g) - 1);
    }

    nmod_poly_reverse(ginv, g, nmod_poly_length(g));
    nmod_poly_inv_series(ginv, ginv, nmod_poly_length(g));

    _nmod_poly_powers_mod_preinv_naive(res_arr, f->coeffs, f->length, n,
                 g->coeffs, g->length, ginv->coeffs, ginv->length, g->mod);

    for (i = 0; i < n; i++)
       _nmod_poly_normalise(res + i);

    nmod_poly_clear(ginv);
    flint_free(res_arr);
}
