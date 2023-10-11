/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2010, 2012 Sebastian Pancratz
    Copyright (C) 2011, 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "gr_poly.h"

void
_nmod_poly_compose(mp_ptr res, mp_srcptr poly1, slong len1,
                               mp_srcptr poly2, slong len2, nmod_t mod)
{
    if (len1 == 1)
        res[0] = poly1[0];
    else if (len2 == 1)
        res[0] = _nmod_poly_evaluate_nmod(poly1, len1, poly2[0], mod);
    else if (len1 <= 7)
        _nmod_poly_compose_horner(res, poly1, len1, poly2, len2, mod);
    else
    {
        /* delegates to Taylor shift, divconquer */
        /* todo: also incorporate the nmod_poly Taylor shift in gr_poly */
        gr_ctx_t ctx;
        _gr_ctx_init_nmod(ctx, &mod);
        GR_MUST_SUCCEED(_gr_poly_compose(res, poly1, len1, poly2, len2, ctx));
    }
}

void nmod_poly_compose(nmod_poly_t res,
                       const nmod_poly_t poly1, const nmod_poly_t poly2)
{
    const slong len1 = poly1->length;
    const slong len2 = poly2->length;

    if (len1 == 0)
    {
        nmod_poly_zero(res);
    }
    else if (len1 == 1 || len2 == 0)
    {
        nmod_poly_fit_length(res, 1);
        res->coeffs[0] = poly1->coeffs[0];
        res->length = (res->coeffs[0] != 0);
    }
    else
    {
        const slong lenr = (len1 - 1) * (len2 - 1) + 1;

        if (res != poly1 && res != poly2)
        {
            nmod_poly_fit_length(res, lenr);
            _nmod_poly_compose(res->coeffs, poly1->coeffs, len1,
                                                   poly2->coeffs, len2, poly1->mod);
        }
        else
        {
            nmod_poly_t t;
            nmod_poly_init2(t, poly1->mod.n, lenr);
            _nmod_poly_compose(t->coeffs, poly1->coeffs, len1,
                                                 poly2->coeffs, len2, poly1->mod);
            nmod_poly_swap(res, t);
            nmod_poly_clear(t);
        }

        res->length = lenr;
        _nmod_poly_normalise(res);
    }
}

void
_nmod_poly_compose_horner(mp_ptr res, mp_srcptr poly1, slong len1,
                                      mp_srcptr poly2, slong len2, nmod_t mod)
{
    if (len1 == 1)
    {
        res[0] = poly1[0];
    }
    else if (len2 == 1)
    {
        res[0] = _nmod_poly_evaluate_nmod(poly1, len1, poly2[0], mod);
    }
    else if (len1 == 2)
    {
        _nmod_vec_scalar_mul_nmod(res, poly2, len2, poly1[1], mod);
        res[0] = n_addmod(res[0], poly1[0], mod.n);
    }
    else
    {
        const slong alloc = (len1 - 1) * (len2 - 1) + 1;
        slong i = len1 - 1, lenr = len2;
        mp_limb_t *t, *t1, *t2;
        TMP_INIT;

        TMP_START;
        t = TMP_ALLOC(alloc * sizeof(mp_limb_t));

        if (len1 % 2 == 0)
        {
            t1 = res;
            t2 = t;
        }
        else
        {
            t1 = t;
            t2 = res;
        }

        /*  Perform the first two steps as one,
            "res = a(m) * poly2 + a(m-1)".      */
        {
            _nmod_vec_scalar_mul_nmod(t1, poly2, len2, poly1[i], mod);
            i--;
            t1[0] = n_addmod(t1[0], poly1[i], mod.n);
        }
        while (i--)
        {
            _nmod_poly_mul(t2, t1, lenr, poly2, len2, mod);
            lenr += len2 - 1;
            MP_PTR_SWAP(t1, t2);
            t1[0] = n_addmod(t1[0], poly1[i], mod.n);
        }

        TMP_END;
    }
}

void nmod_poly_compose_horner(nmod_poly_t res,
                              const nmod_poly_t poly1, const nmod_poly_t poly2)
{
    const slong len1 = poly1->length;
    const slong len2 = poly2->length;

    if (len1 == 0)
    {
        nmod_poly_zero(res);
    }
    else if (len1 == 1 || len2 == 0)
    {
        nmod_poly_fit_length(res, 1);
        res->coeffs[0] = poly1->coeffs[0];
        res->length = (res->coeffs[0] != 0);
    }
    else
    {
        const slong lenr = (len1 - 1) * (len2 - 1) + 1;

        if (res != poly1 && res != poly2)
        {
            nmod_poly_fit_length(res, lenr);
            _nmod_poly_compose_horner(res->coeffs, poly1->coeffs, len1,
                                                   poly2->coeffs, len2, poly1->mod);
        }
        else
        {
            nmod_poly_t t;
            nmod_poly_init2(t, poly1->mod.n, lenr);
            _nmod_poly_compose_horner(t->coeffs, poly1->coeffs, len1,
                                                 poly2->coeffs, len2, poly1->mod);
            nmod_poly_swap(res, t);
            nmod_poly_clear(t);
        }

        res->length = lenr;
        _nmod_poly_normalise(res);
    }
}

void
_nmod_poly_compose_series(mp_ptr res, mp_srcptr poly1, slong len1,
                            mp_srcptr poly2, slong len2, slong n, nmod_t mod)
{
    gr_ctx_t ctx;
    _gr_ctx_init_nmod(ctx, &mod);
    GR_MUST_SUCCEED(_gr_poly_compose_series(res, poly1, len1, poly2, len2, n, ctx));
}

void
nmod_poly_compose_series(nmod_poly_t res,
                    const nmod_poly_t poly1, const nmod_poly_t poly2, slong n)
{
    slong len1 = poly1->length;
    slong len2 = poly2->length;
    slong lenr;

    if (len2 != 0 && poly2->coeffs[0] != 0)
    {
        flint_printf("Exception (nmod_poly_compose_series). Inner polynomial \n"
               "must have zero constant term.\n");
        flint_abort();
    }

    if (len1 == 0 || n == 0)
    {
        nmod_poly_zero(res);
        return;
    }

    if (len2 == 0 || len1 == 1)
    {
        nmod_poly_fit_length(res, 1);
        res->coeffs[0] = poly1->coeffs[0];
        res->length = 1;
        _nmod_poly_normalise(res);
        return;
    }

    lenr = FLINT_MIN((len1 - 1) * (len2 - 1) + 1, n);
    len1 = FLINT_MIN(len1, lenr);
    len2 = FLINT_MIN(len2, lenr);

    if ((res != poly1) && (res != poly2))
    {
        nmod_poly_fit_length(res, lenr);
        _nmod_poly_compose_series(res->coeffs, poly1->coeffs, len1,
                                        poly2->coeffs, len2, lenr, res->mod);
        res->length = lenr;
        _nmod_poly_normalise(res);
    }
    else
    {
        nmod_poly_t t;
        nmod_poly_init2_preinv(t, res->mod.n, res->mod.ninv, lenr);
        _nmod_poly_compose_series(t->coeffs, poly1->coeffs, len1,
                                        poly2->coeffs, len2, lenr, res->mod);
        t->length = lenr;
        _nmod_poly_normalise(t);
        nmod_poly_swap(res, t);
        nmod_poly_clear(t);
    }
}
