/*
    Copyright (C) 2026 Vincent Neiger, Kevin Tran

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "nmod_poly.h"

void nmod_poly_extrapolate_geometric_precomp(nn_ptr oval, slong olen,
                                             nn_srcptr ival, slong ilen,
                                             slong offset,
                                             const nmod_geometric_progression_t G)
{
    /* precomputation has been done */
    FLINT_ASSERT((G->function & UWORD(4)) == UWORD(4));
    /* input/output points are disjoint, and stay within precomputed data length */
    FLINT_ASSERT((offset >= ilen && G->len >= offset+olen)
                 || (offset <= -olen && G->len >= ilen-offset));

    if (olen == 0)
        return;

    if (ilen == 0)
    {
        for (slong i = 0; i < olen; i++)
            oval[i] = 0;
        return;
    }

    if (ilen == 1)
    {
        for (slong i = 0; i < olen; i++)
            oval[i] = ival[0];
        return;
    }

    /* forward extrapolation */
    if (offset > 0)
    {
        TMP_INIT;
        TMP_START;
        nn_ptr tmp = TMP_ALLOC(ilen * sizeof(ulong));

        /* first scaling */
        for (slong i = 0; i < ilen; i++)
        {
            tmp[i] = nmod_mul(G->ext_s3[ilen-1-i], ival[i], G->mod);
            tmp[i] = nmod_mul(G->ext_s2[i], tmp[i], G->mod);
        }

        /* middle product */
        _nmod_poly_mulmid(oval, G->ext_ff->coeffs + offset-ilen, ilen+olen-1, tmp, ilen, ilen-1, ilen+olen-1, G->mod);

        /* second scaling */
        for (slong j = 0; j < olen; j++)
        {
            oval[j] = nmod_mul(G->ext_s2[offset-ilen+j], oval[j], G->mod);
            oval[j] = nmod_mul(G->ext_s1f[offset+j], oval[j], G->mod);
        }
        
        TMP_END;
    }

    /* backward extrapolation */
    else
    {
        TMP_INIT;
        TMP_START;
        nn_ptr tmp = TMP_ALLOC(FLINT_MAX(ilen, olen) * sizeof(ulong));

        /* first scaling */
        for (slong i = 0; i < ilen; i++)
        {
            tmp[i] = nmod_mul(G->ext_s2[ilen-1-i], ival[ilen-1-i], G->mod);
            tmp[i] = nmod_mul(G->ext_s3[i], tmp[i], G->mod);
        }

        /* middle product */
        _nmod_poly_mulmid(oval, G->ext_fb->coeffs - (offset+olen), ilen+olen-1, tmp, ilen, ilen-1, ilen+olen-1, G->mod);

        /* second scaling */
        for (slong j = 0; j < olen; j++)
            tmp[j] = oval[olen - 1 - j];
        for (slong j = 0; j < olen; j++)
        {
            oval[j] = nmod_mul(G->ext_s3[-offset-1-j], tmp[j], G->mod);
            oval[j] = nmod_mul(G->ext_s1b[ilen-1-offset-j], oval[j], G->mod);
        }

        TMP_END;
    }
}

void nmod_poly_extrapolate_geometric(nn_ptr oval, slong olen,
                                     nn_srcptr ival, slong ilen,
                                     slong offset, ulong r, nmod_t mod)
{
    if (olen == 0)
        return;

    if (ilen == 0)
    {
        for (slong i = 0; i < olen; i++)
            oval[i] = 0;
        return;
    }

    if (ilen == 1)
    {
        for (slong i = 0; i < olen; i++)
            oval[i] = ival[0];
        return;
    }

    nmod_geometric_progression_t G;
    slong len = (offset > 0) ? offset+olen : ilen-offset;
    _nmod_geometric_progression_init_function(G, r, len, mod, UWORD(4));
    nmod_poly_extrapolate_geometric_precomp(oval, olen, ival, ilen, offset, G);
    nmod_geometric_progression_clear(G);
}
