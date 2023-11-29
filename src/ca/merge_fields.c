/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca.h"
#include "ca_ext.h"
#include "ca_field.h"
#include "fmpz_mpoly.h"

void
_nf_elem_get_fmpz_poly_den_shallow(fmpz_poly_t pol, fmpz_t den, const nf_elem_t a, const nf_t nf)
{
    if (nf->flag & NF_LINEAR)
    {
        pol->coeffs = (fmpz *) LNF_ELEM_NUMREF(a);
        *den = *LNF_ELEM_DENREF(a);
        pol->length = 1;
        if (fmpz_is_zero(pol->coeffs))
            pol->length--;
    }
    else if (nf->flag & NF_QUADRATIC)
    {
        pol->coeffs = (fmpz *) QNF_ELEM_NUMREF(a);
        *den = *QNF_ELEM_DENREF(a);
        pol->length = 2;
        if (fmpz_is_zero(pol->coeffs + 1))
        {
            pol->length--;
            if (fmpz_is_zero(pol->coeffs))
                pol->length--;
        }
    }
    else
    {
        pol->coeffs = (fmpz *) NF_ELEM_NUMREF(a);
        pol->length = NF_ELEM(a)->length;
        *den = *NF_ELEM_DENREF(a);
    }

    pol->alloc = pol->length;
}


void
ca_merge_fields(ca_t resx, ca_t resy, const ca_t x, const ca_t y, ca_ctx_t ctx)
{
    ca_field_srcptr xfield, yfield, field;
    ca_ext_struct ** ext;
    slong *xgen_map, *ygen_map;
    slong xlen, ylen, ext_len;
    slong ext_alloc;
    slong ix, iy;
    int cmp;

    if (CA_IS_SPECIAL(x) || CA_IS_SPECIAL(y))
    {
        flint_throw(FLINT_ERROR, "ca_merge_fields: inputs must be field elements, not special values\n");
    }

    xfield = CA_FIELD(x, ctx);
    yfield = CA_FIELD(y, ctx);

    if (xfield == yfield || CA_FIELD_IS_QQ(xfield) || CA_FIELD_IS_QQ(yfield))
    {
        ca_set(resx, x, ctx);
        ca_set(resy, y, ctx);
        return;
    }

    if (x == resx || y == resy)
    {
        flint_throw(FLINT_ERROR, "ca_merge_fields: aliasing not implemented!\n");
    }

    xlen = CA_FIELD_LENGTH(xfield);
    ylen = CA_FIELD_LENGTH(yfield);

    ext_alloc = xlen + ylen;
    ext = flint_malloc(ext_alloc * sizeof(ca_ext_struct *));
    ext_len = 0;

    xgen_map = flint_malloc(xlen * sizeof(slong));
    ygen_map = flint_malloc(ylen * sizeof(slong));

/*
    printf("merge fields of len %ld and len %ld\n", xlen, ylen);
    for (ix = 0; ix < xlen; ix++)
    {
        printf("x: "); ca_ext_print(CA_FIELD_EXT_ELEM(xfield, ix), ctx); printf("\n");
    }
    for (iy = 0; iy < ylen; iy++)
    {
        printf("y: "); ca_ext_print(CA_FIELD_EXT_ELEM(yfield, iy), ctx); printf("\n");
    }
*/

    /* merge field lists */

    ix = iy = 0;
    while (ix < xlen || iy < ylen)
    {
        if (ix < xlen && iy < ylen)
        {
            cmp = ca_ext_cmp_repr(CA_FIELD_EXT_ELEM(xfield, ix), CA_FIELD_EXT_ELEM(yfield, iy), ctx);
            cmp = -cmp;  /* more complex first, for elimination order */

            if (cmp == 0)
            {
                if (CA_FIELD_EXT_ELEM(xfield, ix) != CA_FIELD_EXT_ELEM(yfield, iy))
                    flint_throw(FLINT_ERROR, "(%s)\n", __func__);

                ext[ext_len] = CA_FIELD_EXT_ELEM(xfield, ix);
                xgen_map[ix] = ext_len;
                ygen_map[iy] = ext_len;
                ix++;
                iy++;
            }
            else if (cmp < 0)
            {
                ext[ext_len] = CA_FIELD_EXT_ELEM(xfield, ix);
                xgen_map[ix] = ext_len;
                ix++;
            }
            else
            {
                ext[ext_len] = CA_FIELD_EXT_ELEM(yfield, iy);
                ygen_map[iy] = ext_len;
                iy++;
            }

            ext_len++;
        }
        else if (ix < xlen)
        {
            ext[ext_len] = CA_FIELD_EXT_ELEM(xfield, ix);
            xgen_map[ix] = ext_len;
            ix++;
            ext_len++;
        }
        else
        {
            ext[ext_len] = CA_FIELD_EXT_ELEM(yfield, iy);
            ygen_map[iy] = ext_len;
            iy++;
            ext_len++;
        }
    }

/*
    printf("merged length %ld\n", ext_len);
*/

    field = ca_field_cache_insert_ext(CA_CTX_FIELD_CACHE(ctx), ext, ext_len, ctx);

/*
    printf("MERGE FIELDS:\n");
    if (CA_FIELD_LENGTH(xfield) > 100) flint_abort();
    ca_field_print(xfield, ctx); printf("\n");
    if (CA_FIELD_LENGTH(yfield) > 100) flint_abort();
    ca_field_print(yfield, ctx); printf("\n");
    if (CA_FIELD_LENGTH(field) > 100) flint_abort();
    ca_field_print(field, ctx); printf("\n\n");
*/

    if (xfield == field)
    {
        ca_set(resx, x, ctx);
    }
    else
    {
        /* todo: allow aliasing */
        _ca_make_field_element(resx, field, ctx);

        if (CA_FIELD_IS_NF(xfield))
        {
            fmpz_poly_t pol;
            fmpz_t den;

            _nf_elem_get_fmpz_poly_den_shallow(pol, den, CA_NF_ELEM(x), CA_FIELD_NF(xfield));

            fmpz_mpoly_set_gen_fmpz_poly(fmpz_mpoly_q_numref(CA_MPOLY_Q(resx)), xgen_map[0], pol, CA_FIELD_MCTX(field, ctx));
            fmpz_mpoly_set_fmpz(fmpz_mpoly_q_denref(CA_MPOLY_Q(resx)), den, CA_FIELD_MCTX(field, ctx));
        }
        else
        {
            fmpz_mpoly_compose_fmpz_mpoly_gen(fmpz_mpoly_q_numref(CA_MPOLY_Q(resx)),
                                              fmpz_mpoly_q_numref(CA_MPOLY_Q(x)),
                                                xgen_map,
                                                CA_FIELD_MCTX(xfield, ctx),
                                                CA_FIELD_MCTX(field, ctx));

            fmpz_mpoly_compose_fmpz_mpoly_gen(fmpz_mpoly_q_denref(CA_MPOLY_Q(resx)),
                                              fmpz_mpoly_q_denref(CA_MPOLY_Q(x)),
                                                xgen_map,
                                                CA_FIELD_MCTX(xfield, ctx),
                                                CA_FIELD_MCTX(field, ctx));
        }
    }

    if (yfield == field)
    {
        ca_set(resy, y, ctx);
    }
    else
    {
        /* todo: allow aliasing */
        _ca_make_field_element(resy, field, ctx);

        if (CA_FIELD_IS_NF(yfield))
        {
            fmpz_poly_t pol;
            fmpz_t den;

            _nf_elem_get_fmpz_poly_den_shallow(pol, den, CA_NF_ELEM(y), CA_FIELD_NF(yfield));

            fmpz_mpoly_set_gen_fmpz_poly(fmpz_mpoly_q_numref(CA_MPOLY_Q(resy)), ygen_map[0], pol, CA_FIELD_MCTX(field, ctx));
            fmpz_mpoly_set_fmpz(fmpz_mpoly_q_denref(CA_MPOLY_Q(resy)), den, CA_FIELD_MCTX(field, ctx));
        }
        else
        {
            fmpz_mpoly_compose_fmpz_mpoly_gen(fmpz_mpoly_q_numref(CA_MPOLY_Q(resy)),
                                              fmpz_mpoly_q_numref(CA_MPOLY_Q(y)),
                                                ygen_map,
                                                CA_FIELD_MCTX(yfield, ctx),
                                                CA_FIELD_MCTX(field, ctx));

            fmpz_mpoly_compose_fmpz_mpoly_gen(fmpz_mpoly_q_denref(CA_MPOLY_Q(resy)),
                                              fmpz_mpoly_q_denref(CA_MPOLY_Q(y)),
                                                ygen_map,
                                                CA_FIELD_MCTX(yfield, ctx),
                                                CA_FIELD_MCTX(field, ctx));
        }
    }

    flint_free(ext);
    flint_free(xgen_map);
    flint_free(ygen_map);
}
