/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca.h"
#include "ca_ext.h"
#include "ca_field.h"


static int
_ext_vec_equal(ca_ext_struct ** a, ca_ext_struct ** b, slong len)
{
    slong i;

    for (i = 0; i < len; i++)
        if (a[i] != b[i])
            return 0;

    return 1;
}

void
fmpz_mpoly_set_gen_fmpz_poly(fmpz_mpoly_t res, slong var, const fmpz_poly_t pol, const fmpz_mpoly_ctx_t ctx);

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
    slong xfield, yfield, field;
    ca_ext_struct ** ext;
    slong *xgen_map, *ygen_map;
    slong xlen, ylen, ext_len;
    slong ext_alloc;
    slong ix, iy, i;
    int cmp;

    xfield = x->field;
    yfield = y->field;

    if (CA_IS_SPECIAL(x) || CA_IS_SPECIAL(y))
    {
        flint_printf("ca_merge_fields: inputs must be field elements, not special values\n");
        flint_abort();
    }

    if (xfield == yfield)
    {
        ca_set(resx, x, ctx);
        ca_set(resy, y, ctx);
        return;
    }

    /* todo: handle rationals */
    /* (will usually be special-cased, but should be supported here for completeness) */
    if (CA_FIELD_LENGTH(ctx->fields + xfield) == 0 ||
        CA_FIELD_LENGTH(ctx->fields + yfield) == 0)
    {
        flint_printf("QQ in merge_fields\n");
        flint_abort();
    }

    xlen = CA_FIELD_LENGTH(ctx->fields + xfield);
    ylen = CA_FIELD_LENGTH(ctx->fields + yfield);

    ext_alloc = xlen + ylen;
    ext = flint_malloc(ext_alloc * sizeof(ca_ext_struct *));
    ext_len = 0;

    xgen_map = flint_malloc(xlen * sizeof(slong));
    ygen_map = flint_malloc(ylen * sizeof(slong));

/*
    printf("merge fields of len %ld and len %ld\n", xlen, ylen);
*/

    /* merge field lists */

    ix = iy = 0;
    while (ix < xlen || iy < ylen)
    {
        if (ix < xlen && iy < ylen)
        {
            cmp = ca_ext_cmp_repr(CA_FIELD_EXT_ELEM(ctx->fields + xfield, ix), CA_FIELD_EXT_ELEM(ctx->fields + yfield, iy), ctx);
            cmp = -cmp;  /* more complex first, for elimination order */

            if (cmp == 0)
            {
                ext[ext_len] = CA_FIELD_EXT_ELEM(ctx->fields + xfield, ix);
                xgen_map[ix] = ext_len;
                ygen_map[iy] = ext_len;
                ix++;
                iy++;
            }
            else if (cmp == -1)
            {
                ext[ext_len] = CA_FIELD_EXT_ELEM(ctx->fields + xfield, ix);
                xgen_map[ix] = ext_len;
                ix++;
            }
            else
            {
                ext[ext_len] = CA_FIELD_EXT_ELEM(ctx->fields + yfield, iy);
                ygen_map[iy] = ext_len;
                iy++;
            }

            ext_len++;
        }
        else if (ix < xlen)
        {
            ext[ext_len] = CA_FIELD_EXT_ELEM(ctx->fields + xfield, ix);
            xgen_map[ix] = ext_len;
            ix++;
            ext_len++;
        }
        else
        {
            ext[ext_len] = CA_FIELD_EXT_ELEM(ctx->fields + yfield, iy);
            ygen_map[iy] = ext_len;
            iy++;
            ext_len++;
        }
    }

/*
    printf("merge table:\n");
    for (i = 0; i < fields_len; i++)
    {
        flint_printf("%ld\n", fields[i]);
    }
    flint_printf("\n");
*/

    /* check if already cached (todo: needs fast search table) */
    field = -1;
    for (i = 0; i < ctx->fields_len; i++)
    {
        if (CA_FIELD_LENGTH(ctx->fields + i) == ext_len)
        {
            if (_ext_vec_equal(ext, CA_FIELD_EXT(ctx->fields + i), ext_len))
            {
                field = i;
                break;
            }
        }
    }

    if (field == -1)
    {
        field = ctx->fields_len;

        if (field >= ctx->fields_alloc)
        {
            ctx->fields = (ca_field_struct *) flint_realloc(ctx->fields, sizeof(ca_field_struct) * 2 * ctx->fields_alloc);
            ctx->fields_alloc = 2 * ctx->fields_alloc;
        }

        ctx->fields_len = field + 1;
        ca_field_init_multi(ctx->fields + field, ext_len, ctx);

        for (i = 0; i < ext_len; i++)
        {
            ca_field_set_ext(ctx->fields + field, i, ext[i], ctx);
        }

        ca_field_build_ideal(ctx->fields + field, ctx);
    }

/*
    printf("found field\n");
    ca_field_print(ctx->fields + field);
    printf("\n");
*/

    if (xfield == field)
    {
        ca_set(resx, x, ctx);
    }
    else
    {
        /* todo: allow aliasing */
        _ca_make_field_element(resx, field, ctx);

        if (CA_FIELD_IS_NF(ctx->fields + xfield))
        {
            fmpz_poly_t pol;
            fmpz_t den;

            _nf_elem_get_fmpz_poly_den_shallow(pol, den, CA_NF_ELEM(x), CA_FIELD_NF(ctx->fields + xfield));

            fmpz_mpoly_set_gen_fmpz_poly(fmpz_mpoly_q_numref(CA_MPOLY_Q(resx)), xgen_map[0], pol, CA_FIELD_MCTX(ctx->fields + field, ctx));
            fmpz_mpoly_set_fmpz(fmpz_mpoly_q_denref(CA_MPOLY_Q(resx)), den, CA_FIELD_MCTX(ctx->fields + field, ctx));
        }
        else
        {
            fmpz_mpoly_compose_fmpz_mpoly_gen(fmpz_mpoly_q_numref(CA_MPOLY_Q(resx)),
                                              fmpz_mpoly_q_numref(CA_MPOLY_Q(x)),
                                                xgen_map,
                                                CA_FIELD_MCTX(ctx->fields + xfield, ctx),
                                                CA_FIELD_MCTX(ctx->fields + field, ctx));

            fmpz_mpoly_compose_fmpz_mpoly_gen(fmpz_mpoly_q_denref(CA_MPOLY_Q(resx)),
                                              fmpz_mpoly_q_denref(CA_MPOLY_Q(x)),
                                                xgen_map,
                                                CA_FIELD_MCTX(ctx->fields + xfield, ctx),
                                                CA_FIELD_MCTX(ctx->fields + field, ctx));

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

        if (CA_FIELD_IS_NF(ctx->fields + yfield))
        {
            fmpz_poly_t pol;
            fmpz_t den;

            _nf_elem_get_fmpz_poly_den_shallow(pol, den, CA_NF_ELEM(y), CA_FIELD_NF(ctx->fields + yfield));

            fmpz_mpoly_set_gen_fmpz_poly(fmpz_mpoly_q_numref(CA_MPOLY_Q(resy)), ygen_map[0], pol, CA_FIELD_MCTX(ctx->fields + field, ctx));
            fmpz_mpoly_set_fmpz(fmpz_mpoly_q_denref(CA_MPOLY_Q(resy)), den, CA_FIELD_MCTX(ctx->fields + field, ctx));
        }
        else
        {
            fmpz_mpoly_compose_fmpz_mpoly_gen(fmpz_mpoly_q_numref(CA_MPOLY_Q(resy)),
                                              fmpz_mpoly_q_numref(CA_MPOLY_Q(y)),
                                                ygen_map,
                                                CA_FIELD_MCTX(ctx->fields + yfield, ctx),
                                                CA_FIELD_MCTX(ctx->fields + field, ctx));

            fmpz_mpoly_compose_fmpz_mpoly_gen(fmpz_mpoly_q_denref(CA_MPOLY_Q(resy)),
                                              fmpz_mpoly_q_denref(CA_MPOLY_Q(y)),
                                                ygen_map,
                                                CA_FIELD_MCTX(ctx->fields + yfield, ctx),
                                                CA_FIELD_MCTX(ctx->fields + field, ctx));
        }
    }

    flint_free(ext);
    flint_free(xgen_map);
    flint_free(ygen_map);
}

