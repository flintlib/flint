/*=============================================================================

    This file is part of Antic.

    Antic is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version. See <http://www.gnu.org/licenses/>.

=============================================================================*/
/******************************************************************************

    Copyright (C) 2018 Tommy Hofmann

******************************************************************************/

#include "nf_elem.h"

void _nf_elem_get_nmod_poly(nmod_poly_t pol, const nf_elem_t a, const nf_t nf)
{
    if (nf_elem_is_zero(a, nf))
    {
        nmod_poly_zero(pol);
        return;
    }
    if (nf->flag & NF_LINEAR)
    {
        {
            nmod_poly_fit_length(pol, 1);
            pol->coeffs[0] = fmpz_fdiv_ui(LNF_ELEM_NUMREF(a), pol->mod.n);
            _nmod_poly_set_length(pol, 1);
            _nmod_poly_normalise(pol);

        }
    } else if (nf->flag & NF_QUADRATIC)
    {
        nmod_poly_fit_length(pol, 3);
        pol->coeffs[0] = fmpz_fdiv_ui(QNF_ELEM_NUMREF(a), pol->mod.n);
        pol->coeffs[1] = fmpz_fdiv_ui(QNF_ELEM_NUMREF(a) + 1, pol->mod.n);
        pol->coeffs[2] = fmpz_fdiv_ui(QNF_ELEM_NUMREF(a) + 2, pol->mod.n);
        _nmod_poly_set_length(pol, 3);
        _nmod_poly_normalise(pol);
    } else
    {
        slong len = NF_ELEM(a)->length;
        slong i;
        nmod_poly_fit_length(pol, len);
        for (i = 0; i < len; i++)
            pol->coeffs[i] = fmpz_fdiv_ui(NF_ELEM(a)->coeffs + i, pol->mod.n);
        _nmod_poly_set_length(pol, len);
        _nmod_poly_normalise(pol);
    }
}

void nf_elem_get_nmod_poly_den(nmod_poly_t pol, const nf_elem_t a, const nf_t nf, int den)
{
    _nf_elem_get_nmod_poly(pol, a, nf);
    if (den)
    {
        if (nf->flag & NF_LINEAR)
            nmod_poly_scalar_mul_nmod(pol, pol, n_invmod(fmpz_fdiv_ui(LNF_ELEM_DENREF(a), pol->mod.n), pol->mod.n));
        else if (nf->flag & NF_QUADRATIC)
            nmod_poly_scalar_mul_nmod(pol, pol, n_invmod(fmpz_fdiv_ui(QNF_ELEM_DENREF(a), pol->mod.n), pol->mod.n));
        else
            nmod_poly_scalar_mul_nmod(pol, pol, n_invmod(fmpz_fdiv_ui(NF_ELEM_DENREF(a), pol->mod.n), pol->mod.n));
    }
}

void nf_elem_get_nmod_poly(nmod_poly_t pol, const nf_elem_t a, const nf_t nf)
{
    nf_elem_get_nmod_poly_den(pol, a, nf, 1);
}
