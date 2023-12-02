/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "ca.h"
#include "ca_ext.h"
#include "ca_field.h"

typedef struct
{
    ca_ext_ptr * ext;
    const char ** ext_vars;
    slong ext_len;
    ulong flags;
    slong digits;
    int print_where;
}
ca_print_info_t;

void _ca_print(calcium_stream_t out, const ca_t x, ca_print_info_t * info, ca_ctx_t ctx);
void _ca_ext_print(calcium_stream_t out, ca_ext_t x, const char * var, ca_print_info_t * info, ca_ctx_t ctx);

void ca_write(calcium_stream_t out, const ca_t x, ca_ctx_t ctx);

void _ca_all_extensions(ca_ext_ptr ** extensions, slong * length, const ca_t x, ca_ctx_t ctx);

/* strings ********************************************************************/

char * ca_get_str(const ca_t x, ca_ctx_t ctx)
{
    calcium_stream_t out;
    calcium_stream_init_str(out);
    ca_write(out, x, ctx);
    return out->s;
}

/* printing *******************************************************************/

void
ca_ctx_print(ca_ctx_t ctx)
{
    slong i;

    flint_printf("Calcium context with %wd cached fields:\n", CA_CTX_FIELD_CACHE(ctx)->length);
    for (i = 0; i < CA_CTX_FIELD_CACHE(ctx)->length; i++)
    {
        flint_printf("%wd   ", i);
        ca_field_print(CA_CTX_FIELD_CACHE(ctx)->items[i], ctx);
        flint_printf("\n");
    }
    flint_printf("\n");
}

void
calcium_write_fmpz(calcium_stream_t out, const fmpz_t x)
{
    calcium_write_free(out, fmpz_get_str(NULL, 10, x));
}

void
qqbar_write_n(calcium_stream_t out, const qqbar_t x, slong n)
{
    acb_t t;
    slong prec;

    n = FLINT_MAX(1, n);
    prec = n * 3.333 + 10;

    acb_init(t);
    qqbar_get_acb(t, x, prec);

    calcium_write_acb(out, t, n, ARB_STR_NO_RADIUS);
    acb_clear(t);
}

void
calcium_write_nf_elem(calcium_stream_t out,
    const nf_elem_t a, const char * var, const nf_t nf)
{
    const fmpz * num;
    const fmpz * den;
    slong len;

    if (nf_elem_is_zero(a, nf))
    {
        calcium_write(out, "0");
        return;
    }

    if (nf->flag & NF_LINEAR)
    {
        den = LNF_ELEM_DENREF(a);
        num = LNF_ELEM_NUMREF(a);
        len = 1 - fmpz_is_zero(num);
    }
    else if (nf->flag & NF_QUADRATIC)
    {
        den = QNF_ELEM_DENREF(a);
        num = QNF_ELEM_NUMREF(a);
        len = 3;
        while (len != 0 && fmpz_is_zero(num + len - 1))
            len--;
    }
    else
    {
        den = NF_ELEM(a)->den;
        num = NF_ELEM(a)->coeffs;
        len = NF_ELEM(a)->length;
    }

    if (fmpz_is_one(den))
    {
        calcium_write_free(out, _fmpz_poly_get_str_pretty(num, len, var));
    }
    else
    {
        calcium_write(out, "(");
        calcium_write_free(out, _fmpz_poly_get_str_pretty(num, len, var));
        calcium_write(out, ")/");
        calcium_write_fmpz(out, den);
    }
}

void
fmpz_mpoly_q_write_pretty(calcium_stream_t out, const fmpz_mpoly_q_t f, const char ** x, const fmpz_mpoly_ctx_t ctx)
{
    if (fmpz_mpoly_is_one(fmpz_mpoly_q_denref(f), ctx))
    {
        calcium_write_free(out, fmpz_mpoly_get_str_pretty(fmpz_mpoly_q_numref(f), x, ctx));
    }
    else if (fmpz_mpoly_is_fmpz(fmpz_mpoly_q_denref(f), ctx))
    {
        calcium_write(out, "(");
        calcium_write_free(out, fmpz_mpoly_get_str_pretty(fmpz_mpoly_q_numref(f), x, ctx));
        calcium_write(out, ")/");
        calcium_write_free(out, fmpz_mpoly_get_str_pretty(fmpz_mpoly_q_denref(f), x, ctx));
    }
    else
    {
        calcium_write(out, "(");
        calcium_write_free(out, fmpz_mpoly_get_str_pretty(fmpz_mpoly_q_numref(f), x, ctx));
        calcium_write(out, ")/(");
        calcium_write_free(out, fmpz_mpoly_get_str_pretty(fmpz_mpoly_q_denref(f), x, ctx));
        calcium_write(out, ")");
    }
}

void
_ca_field_print(calcium_stream_t out, const ca_field_t K, ca_print_info_t * info, ca_ctx_t ctx)
{
    slong i, j, len, ideal_len;
    const char ** field_vars;

    calcium_write(out, "QQ");

    len = CA_FIELD_LENGTH(K);

    if (len == 0)
        return;

    field_vars = flint_malloc(sizeof(char *) * len);
    for (i = 0; i < len; i++)
        field_vars[i] = "<UNNAMED VARIABLE>";

    j = 0;

    for (i = 0; i < len; i++)
    {
        for ( ; j < info->ext_len; j++)
        {
            if (info->ext[j] == CA_FIELD_EXT_ELEM(K, i))
            {
                field_vars[i] = info->ext_vars[j];
                break;
            }
        }

        if (j == info->ext_len)
        {
            flint_throw(FLINT_ERROR, "_ca_field_print: ext not found!\n");
        }
    }

    calcium_write(out, "(");
    for (i = 0; i < len; i++)
    {
        calcium_write(out, field_vars[i]);
        if (i < len - 1)
            calcium_write(out, ",");
    }
    calcium_write(out, ")");

    if (CA_FIELD_IS_NF(K))
    {
        calcium_write(out, "/<");
        calcium_write_free(out, fmpz_poly_get_str_pretty(QQBAR_POLY(CA_EXT_QQBAR(CA_FIELD_EXT_ELEM(K, 0))), field_vars[0]));
        calcium_write(out, ">");
    }
    else
    {
        ideal_len = CA_FIELD_IDEAL_LENGTH(K);

        if (ideal_len > 0)
        {
            const fmpz_mpoly_ctx_struct * mctx;

            mctx = CA_FIELD_MCTX(K, ctx);

            calcium_write(out, "/<");
            for (i = 0; i < ideal_len; i++)
            {
                calcium_write_free(out, fmpz_mpoly_get_str_pretty(CA_FIELD_IDEAL_ELEM(K, i), (const char **) field_vars, mctx));
                if (i < ideal_len - 1)
                    calcium_write(out, ", ");
            }
            calcium_write(out, ">");
        }
    }

    flint_free(field_vars);
}

void
_ca_ext_print(calcium_stream_t out, ca_ext_t x, const char * var, ca_print_info_t * info, ca_ctx_t ctx)
{
    if (x->head == CA_QQBar)
    {
        if (info->flags & CA_PRINT_N)
        {
            if (qqbar_is_i(CA_EXT_QQBAR(x)))
            {
                calcium_write(out, "I");
            }
            else
            {
                qqbar_write_n(out, CA_EXT_QQBAR(x), info->digits);
            }

            calcium_write(out, " ");
        }

        calcium_write(out, "[");
        calcium_write_free(out, fmpz_poly_get_str_pretty(QQBAR_POLY(CA_EXT_QQBAR(x)), var));
        calcium_write(out, "=0]");
    }
    else
    {
        acb_t t;

        if (info->flags & CA_PRINT_N)
        {
            acb_init(t);
            ca_ext_get_acb_raw(t, x, info->digits * 3.33 + 64, ctx);
            calcium_write_acb(out, t, info->digits, ARB_STR_NO_RADIUS);
            acb_clear(t);
        }

        if (info->flags & CA_PRINT_N)
            calcium_write(out, " [");

        calcium_write(out, calcium_func_name(CA_EXT_HEAD(x)));

        if (CA_EXT_FUNC_NARGS(x) != 0)
        {
            slong i;
            calcium_write(out, "(");
            for (i = 0; i < CA_EXT_FUNC_NARGS(x); i++)
            {
                _ca_print(out, CA_EXT_FUNC_ARGS(x) + i, info, ctx);
                if (i < CA_EXT_FUNC_NARGS(x) - 1)
                    calcium_write(out, ", ");
            }
            calcium_write(out, ")");
        }

        if (info->flags & CA_PRINT_N)
            calcium_write(out, "]");
    }
}

void
_ca_print(calcium_stream_t out, const ca_t x, ca_print_info_t * info, ca_ctx_t ctx)
{
    ca_field_srcptr xfield;
    int print_where;

    if (CA_IS_SPECIAL(x))
    {
        if (CA_IS_UNDEFINED(x))
        {
            calcium_write(out, "Undefined");
        }
        else if (CA_IS_UNKNOWN(x))
        {
            calcium_write(out, "Unknown");
        }
        else if (CA_IS_UNSIGNED_INF(x))
        {
            calcium_write(out, "UnsignedInfinity");
        }
        else  if (CA_IS_SIGNED_INF(x))
        {
            ca_t sgn;
            sgn->field = x->field & ~CA_SPECIAL;
            sgn->elem = x->elem;

            if (CA_IS_QQ(sgn, ctx))
            {
                if (fmpq_sgn(CA_FMPQ(sgn)) > 0)
                    calcium_write(out, "+Infinity");
                else
                    calcium_write(out, "-Infinity");
            }
            else if (CA_IS_QQ_I(sgn, ctx))
            {
                if (fmpz_sgn(QNF_ELEM_NUMREF(CA_NF_ELEM(sgn)) + 1) > 0)
                    calcium_write(out, "+I * Infinity");
                else
                    calcium_write(out, "-I * Infinity");
            }
            else
            {
                calcium_write(out, "Infinity * (");
                _ca_print(out, sgn, info, ctx);
                calcium_write(out, ")");
            }
        }

        return;
    }

    print_where = info->print_where;
    info->print_where = 0;

    xfield = CA_FIELD(x, ctx);

    if (CA_FIELD_IS_QQ(xfield) && fmpz_is_one(CA_FMPQ_DENREF(x)) && fmpz_cmp_si(CA_FMPQ_NUMREF(x), -999999) >= 0 && fmpz_cmp_si(CA_FMPQ_NUMREF(x), 999999) <= 0)
    {
        calcium_write_free(out, fmpq_get_str(NULL, 10, CA_FMPQ(x)));
        return;
    }

    if (info->flags & CA_PRINT_N)
    {
        acb_t t;
        acb_init(t);
        ca_get_acb_raw(t, x, info->digits * 3.33 + 64, ctx);
        calcium_write_acb(out, t, info->digits, ARB_STR_NO_RADIUS);
        acb_clear(t);
    }

    if (info->flags & CA_PRINT_REPR)
    {
        if (info->flags & CA_PRINT_N)
            calcium_write(out, " {");

        if (CA_FIELD_IS_QQ(xfield))
        {
            calcium_write_free(out, fmpq_get_str(NULL, 10, CA_FMPQ(x)));
        }
        else if (CA_FIELD_IS_NF(xfield))
        {
            slong i;
            const char * var = NULL;

            for (i = 0; i < info->ext_len; i++)
            {
                if (info->ext[i] == CA_FIELD_EXT_ELEM(xfield, 0))
                {
                    var = info->ext_vars[i];
                    break;
                }
            }

            calcium_write_nf_elem(out, CA_NF_ELEM(x), var, CA_FIELD_NF(xfield));
        }
        else
        {
            const char ** field_vars;
            slong i, j, len;

            len = CA_FIELD_LENGTH(xfield);
            field_vars = flint_malloc(sizeof(char *) * len);
            for (i = 0; i < len; i++)
                field_vars[i] = "<UNNAMED VARIABLE>";

            j = 0;
            for (i = 0; i < len; i++)
            {
                for ( ; j < info->ext_len; j++)
                {
                    if (info->ext[j] == CA_FIELD_EXT_ELEM(xfield, i))
                    {
                        field_vars[i] = info->ext_vars[j];
                        break;
                    }
                }

                if (j == info->ext_len)
                {
                    flint_throw(FLINT_ERROR, "_ca_field_print: ext not found!\n");
                }
            }

            fmpz_mpoly_q_write_pretty(out, CA_MPOLY_Q(x), (const char **) field_vars, CA_FIELD_MCTX(xfield, ctx));

            flint_free(field_vars);
        }

        if (info->flags & CA_PRINT_FIELD)
        {
            calcium_write(out, "  in  ");
            _ca_field_print(out, xfield, info, ctx);
        }

        if (print_where && info->ext_len > 0)
        {
            slong i, len;
            len = info->ext_len;
            calcium_write(out, " where ");
            for (i = 0; i < len; i++)
            {
                calcium_write(out, info->ext_vars[i]);
                calcium_write(out, " = ");
                _ca_ext_print(out, info->ext[i], info->ext_vars[i], info, ctx);
                if (i != len - 1)
                    calcium_write(out, ", ");
            }
        }

        if (info->flags & CA_PRINT_N)
            calcium_write(out, "}");
    }

/* todo: move... */
}


/* todo: something that doesn't run in quadratic time */

void
_ca_ext_insert_extension(ca_ext_ptr ** extensions, slong * length, ca_ext_t x, ca_ctx_t ctx)
{
    slong i, j;

    for (i = 0; i < *length; i++)
    {
        if (extensions[0][i] == x)
            return;
    }

    if (*length == 0)
    {
        extensions[0] = flint_malloc(4 * sizeof(ca_ext_ptr));
        extensions[0][0] = x;
        *length = 1;
    }
    else
    {
        if (((*length + 1) & (*length)) == 0)
            extensions[0] = flint_realloc(extensions[0], ((*length + 1) * 2) * sizeof(ca_ext_ptr));

        for (i = 0; i < *length; i++)
        {
            if (ca_ext_cmp_repr(extensions[0][i], x, ctx) < 0)
            {
                for (j = *length - 1; j >= i; j--)
                    extensions[0][j + 1] = extensions[0][j];

                extensions[0][i] = x;
                break;
            }
            else if (i == *length - 1)
            {
                extensions[0][*length] = x;
            }
        }

        *length = *length + 1;
    }
}

void
_ca_ext_all_extensions(ca_ext_ptr ** extensions, slong * length, ca_ext_t x, ca_ctx_t ctx)
{
    slong i;

    _ca_ext_insert_extension(extensions, length, x, ctx);

    if (!CA_EXT_IS_QQBAR(x))
    {
        for (i = 0; i < CA_EXT_FUNC_NARGS(x); i++)
            _ca_all_extensions(extensions, length, CA_EXT_FUNC_ARGS(x) + i, ctx);
    }
}

void
_ca_all_extensions(ca_ext_ptr ** extensions, slong * length, const ca_t x, ca_ctx_t ctx)
{
    ca_field_srcptr K;
    slong i;

    if (CA_IS_SPECIAL(x))
    {
        if (CA_IS_SIGNED_INF(x))
        {
            ca_t sgn;
            sgn->field = x->field & ~CA_SPECIAL;
            sgn->elem = x->elem;
            _ca_all_extensions(extensions, length, sgn, ctx);
        }

        return;
    }

    K = CA_FIELD(x, ctx);

    for (i = 0; i < CA_FIELD_LENGTH(K); i++)
    {
        _ca_ext_all_extensions(extensions, length, CA_FIELD_EXT_ELEM(K, i), ctx);
    }
}

void
ca_all_extensions(ca_ext_ptr ** extensions, slong * length, const ca_t x, ca_ctx_t ctx)
{
    *extensions = NULL;
    *length = 0;

    _ca_all_extensions(extensions, length, x, ctx);
}

void
ca_write(calcium_stream_t out, const ca_t x, ca_ctx_t ctx)
{
    ca_print_info_t info;

    ca_ext_ptr * ext;
    slong i, len;
    char * all_vars;
    const char ** vars;

    ca_all_extensions(&ext, &len, x, ctx);

    all_vars = flint_malloc(15 * len);
    vars = flint_malloc(sizeof(char *) * len);
    for (i = 0; i < len; i++)
    {
        if (i < 26)
        {
            all_vars[15 * i] = 'a' + i;
            all_vars[15 * i + 1] = '\0';
        }
        else
        {
            all_vars[15 * i] = 'a' + (i % 26);
            flint_sprintf(&all_vars[15 * i] + 1, "%wd", i / 26);
        }

        vars[i] = all_vars + 15 * i;
    }

    info.ext = ext;
    info.ext_len = len;
    info.ext_vars = vars;
    info.flags = ctx->options[CA_OPT_PRINT_FLAGS];
    info.digits = ctx->options[CA_OPT_PRINT_FLAGS] / CA_PRINT_DIGITS;
    if (info.digits == 0)
        info.digits = 6;
    info.print_where = 1;

    _ca_print(out, x, &info, ctx);

    flint_free(all_vars);
    flint_free(vars);
    flint_free(ext);
}

void
ca_print(const ca_t x, ca_ctx_t ctx)
{
#if 1
    char * s = ca_get_str(x, ctx);
    flint_printf("%s", s);
    flint_free(s);
#else
    calcium_stream_t out;
    calcium_stream_init_file(out, stdout);
    ca_write(out, x, ctx);
#endif
}

void
ca_fprint(FILE * fp, const ca_t x, ca_ctx_t ctx)
{
    calcium_stream_t out;
    calcium_stream_init_file(out, fp);
    ca_write(out, x, ctx);
}

void
ca_printn(const ca_t x, slong n, ca_ctx_t ctx)
{
    ulong save_flags;
    save_flags = ctx->options[CA_OPT_PRINT_FLAGS];
    ctx->options[CA_OPT_PRINT_FLAGS] = CA_PRINT_N | (CA_PRINT_DIGITS * n);
    ca_print(x, ctx);
    ctx->options[CA_OPT_PRINT_FLAGS] = save_flags;
}

void
ca_factor_print(const ca_factor_t fac, ca_ctx_t ctx)
{
    slong i;

    for (i = 0; i < fac->length; i++)
    {
        flint_printf("(");
        ca_print(fac->base + i, ctx); flint_printf(") ^ (");
        ca_print(fac->exp + i, ctx); flint_printf(")\n");
    }
}
