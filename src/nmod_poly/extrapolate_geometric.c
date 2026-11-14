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
#include "nmod_poly/impl.h"
#include "nmod_vec.h"

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

        /* first scaling: tmp[i] = s3[ilen-1-i] ival[i] s2[i] */
        _nmod_poly_reverse(tmp, G->ext_s3, ilen, ilen);
        _nmod_vec_mul(tmp, tmp, ival, ilen, G->mod);
        _nmod_vec_mul(tmp, tmp, G->ext_s2, ilen, G->mod);

        /* middle product; with the stored transform of the whole ext_ff,  */
        /* the coefficients are [offset-1, offset-1+olen) of tmp * ext_ff   */
        int done = (G->ext_ff_fft != NULL)
            && _nmod_geometric_progression_mulmid_fft(oval, offset-1, offset-1+olen, tmp, ilen,
                                                      G->ext_ff_fft, offset-ilen, ilen+olen-1, G->mod);
        if (!done)
            _nmod_poly_mulmid(oval, G->ext_ff->coeffs + offset-ilen, ilen+olen-1, tmp, ilen, ilen-1, ilen+olen-1, G->mod);

        /* second scaling */
        _nmod_vec_mul(oval, oval, G->ext_s2 + offset - ilen, olen, G->mod);
        _nmod_vec_mul(oval, oval, G->ext_s1f + offset, olen, G->mod);

        TMP_END;
    }

    /* backward extrapolation */
    else
    {
        TMP_INIT;
        TMP_START;
        nn_ptr tmp = TMP_ALLOC(FLINT_MAX(ilen, olen) * sizeof(ulong));

        /* first scaling: tmp[i] = s2[ilen-1-i] ival[ilen-1-i] s3[i] */
        _nmod_vec_mul(tmp, G->ext_s2, ival, ilen, G->mod);
        _nmod_poly_reverse(tmp, tmp, ilen, ilen);
        _nmod_vec_mul(tmp, tmp, G->ext_s3, ilen, G->mod);

        /* middle product; with the stored transform of the whole ext_fb,  */
        /* the coefficients are [ilen-1-offset-olen, ilen-1-offset) of tmp * ext_fb */
        int done = (G->ext_fb_fft != NULL)
            && _nmod_geometric_progression_mulmid_fft(oval, ilen-1-offset-olen, ilen-1-offset, tmp, ilen,
                                                      G->ext_fb_fft, -(offset+olen), ilen+olen-1, G->mod);
        if (!done)
            _nmod_poly_mulmid(oval, G->ext_fb->coeffs - (offset+olen), ilen+olen-1, tmp, ilen, ilen-1, ilen+olen-1, G->mod);

        /* second scaling, then reversal:
           oval[olen-1-k] = s3[-offset-olen+k] s1b[ilen-offset-olen+k] oval[k] */
        _nmod_vec_mul(tmp, oval, G->ext_s3 - offset - olen, olen, G->mod);
        _nmod_vec_mul(tmp, tmp, G->ext_s1b + ilen - offset - olen, olen, G->mod);
        _nmod_poly_reverse(oval, tmp, olen, olen);

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
