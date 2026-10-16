/*
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_mod.h"
#include "mpn_mod/impl.h"

int
_mpn_mod_poly_mulmid_karatsuba(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, slong nlo, slong nhi, slong cutoff, gr_ctx_t ctx)
{
    nn_ptr t;
    slong i, l;
    slong nlimbs, slimbs;
    int norm;
    TMP_INIT;
    TMP_START;

    norm = MPN_MOD_CTX_NORM(ctx);

    len1 = FLINT_MIN(len1, nhi);
    len2 = FLINT_MIN(len2, nhi);

    nlimbs = MPN_MOD_CTX_NLIMBS(ctx);

    if (nlo != 0)
    {
        slong nlo2 = (len1 + len2 - 1) - nlo;

        if (len1 > nlo2)
        {
            slong trunc = len1 - nlo2;
            poly1 += trunc * nlimbs;
            len1 -= trunc;
            nlo -= trunc;
            nhi -= trunc;
        }

        if (len2 > nlo2)
        {
            slong trunc = len2 - nlo2;
            poly2 += trunc * nlimbs;
            len2 -= trunc;
            nlo -= trunc;
            nhi -= trunc;
        }
    }

    if (cutoff == -1)
        cutoff = _flint_mpn_poly_mulmid_cutoff(nlimbs, nlimbs, poly1 == poly2 && len1 == len2);

    slimbs = _flint_mpn_poly_mulmid_slimbs(FLINT_MIN(len1, len2), nlimbs, norm, nlimbs, norm, FLINT_MPN_POLY_MUL_UNSIGNED, cutoff);

    t = TMP_ALLOC(sizeof(ulong) * slimbs * (nhi - nlo));

    _flint_mpn_poly_mulmid_karatsuba(t, poly1, len1, nlimbs, norm, poly2, len2, nlimbs, norm, nlo, nhi, slimbs, cutoff, 0);

    for (i = nlo; i < nhi; i++)
    {
        l = slimbs;
        MPN_NORM(t + (i - nlo) * slimbs, l);
        mpn_mod_set_mpn(res + (i - nlo) * nlimbs, t + (i - nlo) * slimbs, l, ctx);
    }

    TMP_END;
    return GR_SUCCESS;
}

/* Classical or Karatsuba multiplication, chosen by the shared mpn code. */
int
_mpn_mod_poly_mulmid_mpn(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, slong nlo, slong nhi, gr_ctx_t ctx)
{
    nn_ptr t;
    slong i, l, nlimbs, slimbs, cutoff, ilo, ihi;
    int norm;
    TMP_INIT;

    nlimbs = MPN_MOD_CTX_NLIMBS(ctx);

    len1 = FLINT_MIN(len1, nhi);
    len2 = FLINT_MIN(len2, nhi);

    if (nlo != 0)
    {
        slong nlo2 = (len1 + len2 - 1) - nlo;

        if (len1 > nlo2)
        {
            slong trunc = len1 - nlo2;
            poly1 += trunc * nlimbs;
            len1 -= trunc;
            nlo -= trunc;
            nhi -= trunc;
        }

        if (len2 > nlo2)
        {
            slong trunc = len2 - nlo2;
            poly2 += trunc * nlimbs;
            len2 -= trunc;
            nlo -= trunc;
            nhi -= trunc;
        }
    }

    if (nhi == 1)
        return mpn_mod_mul(res, poly1, poly2, ctx);

    if (len1 == 1)
        return _mpn_mod_vec_mul_scalar(res, poly2 + nlo * nlimbs, nhi - nlo, poly1, ctx);

    if (len2 == 1)
        return _mpn_mod_vec_mul_scalar(res, poly1 + nlo * nlimbs, nhi - nlo, poly2, ctx);

    /* the residues are unsigned integers of the size of the modulus, so
       the packing is fixed and only the output size needs to be chosen
       (for the lengths after the truncation above) */
    norm = MPN_MOD_CTX_NORM(ctx);
    cutoff = _flint_mpn_poly_mulmid_cutoff(nlimbs, nlimbs, poly1 == poly2 && len1 == len2);
    slimbs = _flint_mpn_poly_mulmid_slimbs(FLINT_MIN(len1, len2), nlimbs, norm, nlimbs, norm, FLINT_MPN_POLY_MUL_UNSIGNED, cutoff);

    ilo = nlo;
    ihi = nhi;

    /* In the classical range, the first and last coefficients (single
       products) are cheaper with the fixed-size modular multiplication
       than as unreduced products followed by the general reduction.
       (Karatsuba keeps the full window, so that a full product is
       recognized as such.) */
    if (!_flint_mpn_poly_mulmid_use_karatsuba(FLINT_MIN(len1, len2), norm, FLINT_MPN_POLY_MUL_UNSIGNED, cutoff))
    {
        if (nlo == 0)
        {
            mpn_mod_mul(res, poly1, poly2, ctx);
            ilo = 1;
        }

        if (nhi == len1 + len2 - 1)
        {
            mpn_mod_mul(res + (nhi - 1 - nlo) * nlimbs, poly1 + (len1 - 1) * nlimbs, poly2 + (len2 - 1) * nlimbs, ctx);
            ihi = nhi - 1;
        }

        if (ilo >= ihi)
            return GR_SUCCESS;
    }

    TMP_START;
    t = TMP_ALLOC(sizeof(ulong) * slimbs * (ihi - ilo));

    _flint_mpn_poly_mulmid(t, poly1, len1, nlimbs, norm, poly2, len2, nlimbs, norm, ilo, ihi, slimbs, FLINT_MPN_POLY_MUL_UNSIGNED);

    for (i = ilo; i < ihi; i++)
    {
        l = slimbs;
        MPN_NORM(t + (i - ilo) * slimbs, l);
        mpn_mod_set_mpn(res + (i - nlo) * nlimbs, t + (i - ilo) * slimbs, l, ctx);
    }

    TMP_END;
    return GR_SUCCESS;
}

int
_mpn_mod_poly_mullow_karatsuba(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, slong len, slong cutoff, gr_ctx_t ctx)
{
    return _mpn_mod_poly_mulmid_karatsuba(res, poly1, len1, poly2, len2, 0, len, cutoff, ctx);
}

