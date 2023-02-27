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

#if __FLINT_RELEASE >= 20700
void _nf_elem_get_fmpz_mod_poly(fmpz_mod_poly_t pol, const nf_elem_t a,
                                       const nf_t nf, const fmpz_mod_ctx_t ctx)
#else
void _nf_elem_get_fmpz_mod_poly(fmpz_mod_poly_t pol, const nf_elem_t a, const nf_t nf)
#endif
{
    if (nf_elem_is_zero(a, nf))
    {
        FMPZ_MOD_POLY_ZERO(pol, ctx);
        
        return;
    }
    if (nf->flag & NF_LINEAR)
    {
        {
            FMPZ_MOD_POLY_FIT_LENGTH(pol, 1, ctx);
        
            FMPZ_MOD(pol->coeffs + 0, LNF_ELEM_NUMREF(a), ctx, &(pol->p));
        
            _fmpz_mod_poly_set_length(pol, 1);
            _fmpz_mod_poly_normalise(pol);

        }
    } else if (nf->flag & NF_QUADRATIC)
    {
        FMPZ_MOD_POLY_FIT_LENGTH(pol, 3, ctx);
        
        FMPZ_MOD(pol->coeffs + 0, QNF_ELEM_NUMREF(a), ctx, &(pol->p));
        FMPZ_MOD(pol->coeffs + 1, QNF_ELEM_NUMREF(a) + 1, ctx, &(pol->p));
        FMPZ_MOD(pol->coeffs + 2, QNF_ELEM_NUMREF(a) + 2, ctx, &(pol->p));
        
        _fmpz_mod_poly_set_length(pol, 3);
        _fmpz_mod_poly_normalise(pol);
    } else
    {
        slong len = NF_ELEM(a)->length;
        slong i;

        FMPZ_MOD_POLY_FIT_LENGTH(pol, len, ctx);
        
        for (i = 0; i < len; i++)
            FMPZ_MOD(pol->coeffs + i, NF_ELEM_NUMREF(a) + i, ctx, &(pol->p));
        
        _fmpz_mod_poly_set_length(pol, len);
        _fmpz_mod_poly_normalise(pol);
    }
}

#if __FLINT_RELEASE >= 20700
void nf_elem_get_fmpz_mod_poly_den(fmpz_mod_poly_t pol, const nf_elem_t a,
                              const nf_t nf, int den, const fmpz_mod_ctx_t ctx)
#else
void nf_elem_get_fmpz_mod_poly_den(fmpz_mod_poly_t pol, const nf_elem_t a, const nf_t nf, int den)
#endif
{
#if __FLINT_RELEASE >= 20700
    _nf_elem_get_fmpz_mod_poly(pol, a, nf, ctx);
#else
    _nf_elem_get_fmpz_mod_poly(pol, a, nf);
#endif
    if (den)
    {
        if (nf->flag & NF_LINEAR)
            FMPZ_MOD_POLY_SCALAR_DIV_FMPZ(pol, pol, LNF_ELEM_DENREF(a), ctx);
        else if (nf->flag & NF_QUADRATIC)
            FMPZ_MOD_POLY_SCALAR_DIV_FMPZ(pol, pol, QNF_ELEM_DENREF(a), ctx);
        else
            FMPZ_MOD_POLY_SCALAR_DIV_FMPZ(pol, pol, NF_ELEM_DENREF(a), ctx);
    }
}

#if __FLINT_RELEASE >= 20700
void nf_elem_get_fmpz_mod_poly(fmpz_mod_poly_t pol, const nf_elem_t a,
                                       const nf_t nf, const fmpz_mod_ctx_t ctx)
{
    nf_elem_get_fmpz_mod_poly_den(pol, a, nf, 1, ctx);
}
#else
void nf_elem_get_fmpz_mod_poly(fmpz_mod_poly_t pol, const nf_elem_t a, const nf_t nf)
{
    nf_elem_get_fmpz_mod_poly_den(pol, a, nf, 1);
}
#endif
