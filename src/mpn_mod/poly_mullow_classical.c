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

/* currently defined in poly_mullow_karatsuba.c */
static void
mpn_mod_set_mpn2(nn_ptr res, nn_srcptr s, slong l, gr_ctx_t ctx)
{
    MPN_NORM(s, l);
    mpn_mod_set_mpn(res, s, l, ctx);
}

int
_mpn_mod_poly_mulmid_classical(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, slong nlo, slong nhi, gr_ctx_t ctx)
{
    slong i;
    slong nlimbs, slimbs;

    if (nhi == 1)
        return mpn_mod_mul(res, poly1, poly2, ctx);

    len1 = FLINT_MIN(len1, nhi);
    len2 = FLINT_MIN(len2, nhi);

    nlimbs = MPN_MOD_CTX_NLIMBS(ctx);

    if (len1 == 1)
        return _mpn_mod_vec_mul_scalar(res, poly2 + nlo * nlimbs, nhi - nlo, poly1, ctx);

    if (len2 == 1)
        return _mpn_mod_vec_mul_scalar(res, poly1 + nlo * nlimbs, nhi - nlo, poly2, ctx);

    if (nlo == 0)
        mpn_mod_mul(res, poly1, poly2, ctx);

    res -= nlo * nlimbs;

    if (nhi == len1 + len2 - 1)
        mpn_mod_mul(res + (len1 + len2 - 2) * nlimbs, poly1 + (len1 - 1) * nlimbs, poly2 + (len2 - 1) * nlimbs, ctx);

    if (nhi <= 1 || nlo >= len1 + len2 - 2)
        return GR_SUCCESS;

    slimbs = _flint_mpn_poly_mulmid_slimbs(FLINT_MIN(len1, len2), nlimbs, MPN_MOD_CTX_NORM(ctx), nlimbs, MPN_MOD_CTX_NORM(ctx), FLINT_MPN_POLY_MUL_UNSIGNED, WORD_MAX);

    /* The first and last coefficients have already been computed above;
       the remaining ones are computed as unreduced dot products of
       fixed limb size and then reduced. */
    {
        slong ilo = FLINT_MAX(nlo, 1);
        slong ihi = FLINT_MIN(nhi, len1 + len2 - 2);
        nn_ptr t;
        TMP_INIT;

        TMP_START;
        t = TMP_ALLOC(sizeof(ulong) * (ihi - ilo) * slimbs);

        _flint_mpn_poly_mulmid_classical(t, poly1, len1, nlimbs, poly2, len2, nlimbs, ilo, ihi, slimbs, 0);

        for (i = ilo; i < ihi; i++)
            mpn_mod_set_mpn2(res + i * nlimbs, t + (i - ilo) * slimbs, slimbs, ctx);

        TMP_END;
    }

    return GR_SUCCESS;
}

int
_mpn_mod_poly_mullow_classical(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, slong len, gr_ctx_t ctx)
{
    return _mpn_mod_poly_mulmid_classical(res, poly1, len1, poly2, len2, 0, len, ctx);
}

