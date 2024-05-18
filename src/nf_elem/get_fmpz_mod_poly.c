/*
    Copyright (C) 2018 Tommy Hofmann

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_poly.h"
#include "nf_elem.h"

void _nf_elem_get_fmpz_mod_poly(fmpz_mod_poly_t pol, const nf_elem_t a,
                                       const nf_t nf, const fmpz_mod_ctx_t ctx)
{
    if (nf_elem_is_zero(a, nf))
    {
        fmpz_mod_poly_zero(pol, ctx);

        return;
    }
    if (nf->flag & NF_LINEAR)
    {
        {
            fmpz_mod_poly_fit_length(pol, 1, ctx);

            fmpz_mod(pol->coeffs + 0, LNF_ELEM_NUMREF(a), ctx->n);

            _fmpz_mod_poly_set_length(pol, 1);
            _fmpz_mod_poly_normalise(pol);

        }
    } else if (nf->flag & NF_QUADRATIC)
    {
        fmpz_mod_poly_fit_length(pol, 3, ctx);

        fmpz_mod(pol->coeffs + 0, QNF_ELEM_NUMREF(a) + 0, ctx->n);
        fmpz_mod(pol->coeffs + 1, QNF_ELEM_NUMREF(a) + 1, ctx->n);
        fmpz_mod(pol->coeffs + 2, QNF_ELEM_NUMREF(a) + 2, ctx->n);

        _fmpz_mod_poly_set_length(pol, 3);
        _fmpz_mod_poly_normalise(pol);
    } else
    {
        slong len = NF_ELEM(a)->length;
        slong i;

        fmpz_mod_poly_fit_length(pol, len, ctx);

        for (i = 0; i < len; i++)
            fmpz_mod(pol->coeffs + i, NF_ELEM_NUMREF(a) + i, ctx->n);

        _fmpz_mod_poly_set_length(pol, len);
        _fmpz_mod_poly_normalise(pol);
    }
}

void nf_elem_get_fmpz_mod_poly_den(fmpz_mod_poly_t pol, const nf_elem_t a,
                              const nf_t nf, int den, const fmpz_mod_ctx_t ctx)
{
    _nf_elem_get_fmpz_mod_poly(pol, a, nf, ctx);
    if (den)
    {
        if (nf->flag & NF_LINEAR)
            fmpz_mod_poly_scalar_div_fmpz(pol, pol, LNF_ELEM_DENREF(a), ctx);
        else if (nf->flag & NF_QUADRATIC)
            fmpz_mod_poly_scalar_div_fmpz(pol, pol, QNF_ELEM_DENREF(a), ctx);
        else
            fmpz_mod_poly_scalar_div_fmpz(pol, pol, NF_ELEM_DENREF(a), ctx);
    }
}

void nf_elem_get_fmpz_mod_poly(fmpz_mod_poly_t pol, const nf_elem_t a,
                                       const nf_t nf, const fmpz_mod_ctx_t ctx)
{
    nf_elem_get_fmpz_mod_poly_den(pol, a, nf, 1, ctx);
}
