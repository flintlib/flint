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

void _ca_field_ideal_insert_clear_mpoly(ca_field_t K, fmpz_mpoly_t poly, fmpz_mpoly_ctx_t mctx, ca_ctx_t ctx);


/* todo: optimize */
truth_t ca_check_equal_neg(const ca_t x, const ca_t y, ca_ctx_t ctx)
{
    ca_t t;
    truth_t res;
    ca_init(t, ctx);
    ca_neg(t, y, ctx);
    res = ca_check_equal(x, t, ctx);
    ca_clear(t, ctx);
    return res;
}

/* set a*x_a + b*x_b + c */
void
fmpz_mpoly_set_linear_three_term_si(fmpz_mpoly_t poly, slong a, slong xa, slong b, slong xb, slong c, const fmpz_mpoly_ctx_t ctx)
{
    ulong * exp;
    exp = flint_calloc(ctx->minfo->nvars, sizeof(ulong));

    if (xa == xb)
    {
        flint_throw(FLINT_ERROR, "fmpz_mpoly_set_linear_three_term_si\n");
    }

    fmpz_mpoly_set_si(poly, c, ctx);

    exp[xa] = 1;
    fmpz_mpoly_set_coeff_si_ui(poly, a, exp, ctx);
    exp[xa] = 0;

    exp[xb] = 1;
    fmpz_mpoly_set_coeff_si_ui(poly, b, exp, ctx);

    flint_free(exp);
}

/* set a*x_a*x_a2 + b*x_b + c */
void
fmpz_mpoly_set_linear2_three_term_si(fmpz_mpoly_t poly, slong a, slong xa, slong xa2, slong b, slong xb, slong c, const fmpz_mpoly_ctx_t ctx)
{
    ulong * exp;
    exp = flint_calloc(ctx->minfo->nvars, sizeof(ulong));

    if (xa == xb || xa == xa2)
    {
        flint_throw(FLINT_ERROR, "fmpz_mpoly_set_linear2_three_term_si\n");
    }

    fmpz_mpoly_set_si(poly, c, ctx);

    exp[xa] = 1;
    exp[xa2] = 1;
    fmpz_mpoly_set_coeff_si_ui(poly, a, exp, ctx);
    exp[xa] = 0;
    exp[xa2] = 0;

    exp[xb] = 1;
    fmpz_mpoly_set_coeff_si_ui(poly, b, exp, ctx);

    flint_free(exp);
}

/* Set the term c * x_var^x_exp */
void
fmpz_mpoly_set_coeff_si_x(fmpz_mpoly_t poly,
        slong c,
        slong x_var, slong x_exp,
        const fmpz_mpoly_ctx_t ctx)
{
    ulong * exp;
    slong i, len;
    TMP_INIT;

    len = ctx->minfo->nvars;

    TMP_START;
    exp = TMP_ALLOC(len * sizeof(ulong));
    for (i = 0; i < len; i++)
        exp[i] = 0;

    exp[x_var] = x_exp;
    fmpz_mpoly_set_coeff_si_ui(poly, c, exp, ctx);
    TMP_END;
}


void
fmpz_mpoly_set_coeff_si_xy(fmpz_mpoly_t poly,
        slong c,
        slong x_var, ulong x_exp,
        slong y_var, ulong y_exp,
        const fmpz_mpoly_ctx_t ctx)
{
    ulong * exp;
    slong i, len;
    TMP_INIT;

    len = ctx->minfo->nvars;

    TMP_START;
    exp = TMP_ALLOC(len * sizeof(ulong));
    for (i = 0; i < len; i++)
        exp[i] = 0;

    exp[x_var] = x_exp;
    exp[y_var] = y_exp;
    fmpz_mpoly_set_coeff_si_ui(poly, c, exp, ctx);
    TMP_END;
}

static void
ideal_mixed_erfi(ca_field_t K, slong i, slong j, int have_i, slong index_i, ca_ctx_t ctx)
{
    calcium_func_code Fi;
    ca_t ix;
    ca_ptr x, y;
    fmpz_mpoly_t poly;
    const fmpz_mpoly_ctx_struct * mctx = CA_FIELD_MCTX(K, ctx);

    x = CA_EXT_FUNC_ARGS(CA_FIELD_EXT_ELEM(K, i));
    y = CA_EXT_FUNC_ARGS(CA_FIELD_EXT_ELEM(K, j));
    Fi = CA_EXT_HEAD(CA_FIELD_EXT_ELEM(K, i));

    ca_init(ix, ctx);
    ca_i(ix, ctx);
    ca_mul(ix, ix, x, ctx);

    if (ca_check_equal(ix, y, ctx) == T_TRUE)
    {
        if (have_i)
        {
            fmpz_mpoly_init(poly, CA_FIELD_MCTX(K, ctx));
            if (Fi == CA_Erf)
            {
                /* erf(x) + i*erfi(i*x) */
                fmpz_mpoly_set_coeff_si_x(poly, 1, i, 1, mctx);
                fmpz_mpoly_set_coeff_si_xy(poly, 1, j, 1, index_i, 1, mctx);
            }
            else
            {
                /* erfc(x) - i*erfi(i*x) - 1 */
                fmpz_mpoly_set_si(poly, -1, mctx);
                fmpz_mpoly_set_coeff_si_x(poly, 1, i, 1, mctx);
                fmpz_mpoly_set_coeff_si_xy(poly, -1, j, 1, index_i, 1, mctx);
            }
            _ca_field_ideal_insert_clear_mpoly(K, poly, CA_FIELD_MCTX(K, ctx), ctx);
        }

        if (have_i)
        {
            fmpz_mpoly_init(poly, CA_FIELD_MCTX(K, ctx));
            if (Fi == CA_Erf)
            {
                /* i*erf(x) - erfi(i*x) */
                fmpz_mpoly_set_coeff_si_xy(poly, 1, i, 1, index_i, 1, mctx);
                fmpz_mpoly_set_coeff_si_x(poly, -1, j, 1, mctx);
            }
            else
            {
                /* i*erfc(x) + erfi(i*x) - i */
                fmpz_mpoly_set_coeff_si_xy(poly, 1, i, 1, index_i, 1, mctx);
                fmpz_mpoly_set_coeff_si_x(poly, 1, j, 1, mctx);
                fmpz_mpoly_set_coeff_si_x(poly, -1, index_i, 1, mctx);
            }
            _ca_field_ideal_insert_clear_mpoly(K, poly, CA_FIELD_MCTX(K, ctx), ctx);
        }

        fmpz_mpoly_init(poly, CA_FIELD_MCTX(K, ctx));
        if (Fi == CA_Erf)
        {
            /* erf(x)^2 + erfi(i*x)^2 = 0 */
            fmpz_mpoly_set_coeff_si_x(poly, 1, i, 2, mctx);
            fmpz_mpoly_set_coeff_si_x(poly, 1, j, 2, mctx);
        }
        else
        {
            /* erfc(x)^2 - 2*erfc(x) + erfi(i*x)**2 + 1 */
            fmpz_mpoly_set_si(poly, 1, mctx);
            fmpz_mpoly_set_coeff_si_x(poly, 1, i, 2, mctx);
            fmpz_mpoly_set_coeff_si_x(poly, -2, i, 1, mctx);
            fmpz_mpoly_set_coeff_si_x(poly, 1, j, 2, mctx);
        }
        _ca_field_ideal_insert_clear_mpoly(K, poly, CA_FIELD_MCTX(K, ctx), ctx);
    }
    else if (ca_check_equal_neg(ix, y, ctx) == T_TRUE)
    {
        if (have_i)
        {
            fmpz_mpoly_init(poly, CA_FIELD_MCTX(K, ctx));
            if (Fi == CA_Erf)
            {
                /* erf(x) - i*erfi(-i*x) */
                fmpz_mpoly_set_coeff_si_x(poly, 1, i, 1, mctx);
                fmpz_mpoly_set_coeff_si_xy(poly, -1, j, 1, index_i, 1, mctx);
            }
            else
            {
                /* erfc(x) + i*erfi(-i*x) - 1 */
                fmpz_mpoly_set_si(poly, -1, mctx);
                fmpz_mpoly_set_coeff_si_x(poly, 1, i, 1, mctx);
                fmpz_mpoly_set_coeff_si_xy(poly, 1, j, 1, index_i, 1, mctx);
            }
            _ca_field_ideal_insert_clear_mpoly(K, poly, CA_FIELD_MCTX(K, ctx), ctx);
        }

        if (have_i)
        {
            fmpz_mpoly_init(poly, CA_FIELD_MCTX(K, ctx));
            if (Fi == CA_Erf)
            {
                /* i*erf(x) + erfi(-i*x) */
                fmpz_mpoly_set_coeff_si_xy(poly, 1, i, 1, index_i, 1, mctx);
                fmpz_mpoly_set_coeff_si_x(poly, 1, j, 1, mctx);
            }
            else
            {
                /* i*erfc(x) - erfi(-i*x) - i */
                fmpz_mpoly_set_coeff_si_xy(poly, 1, i, 1, index_i, 1, mctx);
                fmpz_mpoly_set_coeff_si_x(poly, -1, j, 1, mctx);
                fmpz_mpoly_set_coeff_si_x(poly, -1, index_i, 1, mctx);
            }
            _ca_field_ideal_insert_clear_mpoly(K, poly, CA_FIELD_MCTX(K, ctx), ctx);
        }

        fmpz_mpoly_init(poly, CA_FIELD_MCTX(K, ctx));
        if (Fi == CA_Erf)
        {
            /* erf(x)^2 + erfi(-i*x)^2 = 0 */
            fmpz_mpoly_set_coeff_si_x(poly, 1, i, 2, mctx);
            fmpz_mpoly_set_coeff_si_x(poly, 1, j, 2, mctx);
        }
        else
        {
            /* erfc(x)^2 - 2*erfc(x) + erfi(-i*x)**2 + 1 */
            fmpz_mpoly_set_si(poly, 1, mctx);
            fmpz_mpoly_set_coeff_si_x(poly, 1, i, 2, mctx);
            fmpz_mpoly_set_coeff_si_x(poly, -2, i, 1, mctx);
            fmpz_mpoly_set_coeff_si_x(poly, 1, j, 2, mctx);
        }
        _ca_field_ideal_insert_clear_mpoly(K, poly, CA_FIELD_MCTX(K, ctx), ctx);
    }

    ca_clear(ix, ctx);
}

/*
erf(x) - erf(x)            = 0
erfc(x) - erfc(x)          = 0
erfi(x) - erfi(x)          = 0
erf(x) + erfc(x) - 1       = 0

erf(x) + erf(-x)           = 0
erfc(x) + erfc(-x) - 2     = 0
erfi(x) + erfi(-x)         = 0
erf(x) - erfc(-x) + 1      = 0
erfc(x) - erf(-x) - 1      = 0

erf(x) + i*erfi(i*x)       = 0
erf(x) - i*erfi(-i*x)      = 0
erfc(x) - i*erfi(i*x) - 1  = 0
erfc(x) + i*erfi(-i*x) - 1 = 0
*/
void
ca_field_build_ideal_erf(ca_field_t K, ca_ctx_t ctx)
{
    slong i, j, len, num_erf, index_i;
    calcium_func_code Fi, Fj;
    int have_i;

    len = CA_FIELD_LENGTH(K);

    if (len < 2)
        return;

    index_i = 0;
    have_i = 0;
    num_erf = 0;
    for (i = 0; i < len; i++)
    {
        Fi = CA_EXT_HEAD(CA_FIELD_EXT_ELEM(K, i));

        if (Fi == CA_Erf || Fi == CA_Erfc || Fi == CA_Erfi)
        {
            num_erf++;
        }
        else if (CA_FIELD_EXT_ELEM(K, i) == CA_FIELD_EXT_ELEM(ctx->field_qq_i, 0))
        {
            have_i = 1;
            index_i = i;
        }
    }

    if (num_erf >= 2)
    {
        for (i = 0; i < len; i++)
        {
            Fi = CA_EXT_HEAD(CA_FIELD_EXT_ELEM(K, i));

            if (Fi == CA_Erf || Fi == CA_Erfc || Fi == CA_Erfi)
            {
                for (j = i + 1; j < len; j++)
                {
                    Fj = CA_EXT_HEAD(CA_FIELD_EXT_ELEM(K, j));

                    if (Fj == CA_Erf || Fj == CA_Erfc || Fj == CA_Erfi)
                    {
                        if ((Fj == CA_Erfi && (Fi == CA_Erf || Fi == CA_Erfc)))
                        {
                            ideal_mixed_erfi(K, i, j, have_i, index_i, ctx);
                            continue;
                        }

                        if ((Fi == CA_Erfi && (Fj == CA_Erf || Fj == CA_Erfc)))
                        {
                            ideal_mixed_erfi(K, j, i, have_i, index_i, ctx);
                            continue;
                        }

                        if (Fi == Fj || (Fi == CA_Erf && Fj == CA_Erfc) || (Fi == CA_Erfc && Fj == CA_Erf))
                        {
                            if (ca_check_equal(CA_EXT_FUNC_ARGS(CA_FIELD_EXT_ELEM(K, i)),
                                               CA_EXT_FUNC_ARGS(CA_FIELD_EXT_ELEM(K, j)), ctx) == T_TRUE)
                            {
                                /*
                                erf(x) - erf(x)            = 0
                                erfc(x) - erfc(x)          = 0
                                erfi(x) - erfi(x)          = 0
                                erf(x) + erfc(x) - 1       = 0
                                */
                                fmpz_mpoly_t poly;
                                fmpz_mpoly_init(poly, CA_FIELD_MCTX(K, ctx));
                                if (Fi == Fj)
                                    fmpz_mpoly_set_linear_three_term_si(poly, 1, i, -1, j, 0, CA_FIELD_MCTX(K, ctx));
                                else
                                    fmpz_mpoly_set_linear_three_term_si(poly, 1, i, 1, j, -1, CA_FIELD_MCTX(K, ctx));
                                _ca_field_ideal_insert_clear_mpoly(K, poly, CA_FIELD_MCTX(K, ctx), ctx);
                            }
                            else if (ca_check_equal_neg(CA_EXT_FUNC_ARGS(CA_FIELD_EXT_ELEM(K, i)),
                                                        CA_EXT_FUNC_ARGS(CA_FIELD_EXT_ELEM(K, j)), ctx) == T_TRUE)
                            {
                                /*
                                erf(x) + erf(-x)           = 0
                                erfc(x) + erfc(-x) - 2     = 0
                                erfi(x) + erfi(-x)         = 0
                                erf(x) - erfc(-x) + 1      = 0
                                erfc(x) - erf(-x) - 1      = 0
                                */
                                fmpz_mpoly_t poly;
                                fmpz_mpoly_init(poly, CA_FIELD_MCTX(K, ctx));
                                if (Fi == Fj)
                                {
                                    if (Fi == CA_Erfc)
                                        fmpz_mpoly_set_linear_three_term_si(poly, 1, i, 1, j, -2, CA_FIELD_MCTX(K, ctx));
                                    else
                                        fmpz_mpoly_set_linear_three_term_si(poly, 1, i, 1, j, 0, CA_FIELD_MCTX(K, ctx));
                                }
                                else
                                {
                                    if (Fi == CA_Erf)
                                        fmpz_mpoly_set_linear_three_term_si(poly, 1, i, -1, j, 1, CA_FIELD_MCTX(K, ctx));
                                    else
                                        fmpz_mpoly_set_linear_three_term_si(poly, 1, i, -1, j, -1, CA_FIELD_MCTX(K, ctx));
                                }
                                _ca_field_ideal_insert_clear_mpoly(K, poly, CA_FIELD_MCTX(K, ctx), ctx);
                            }
                        }
                    }
                }
            }
        }
    }
}
