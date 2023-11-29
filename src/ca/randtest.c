/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca.h"

void
ca_randtest_same_nf(ca_t res, flint_rand_t state, const ca_t x, slong bits, slong den_bits, ca_ctx_t ctx)
{
    if (n_randint(state, 2) || CA_IS_QQ(x, ctx))
    {
        _ca_make_fmpq(res, ctx);
        fmpz_randtest(CA_FMPQ_NUMREF(res), state, bits);
        fmpz_randtest_not_zero(CA_FMPQ_DENREF(res), state, den_bits);
        fmpz_abs(CA_FMPQ_DENREF(res), CA_FMPQ_DENREF(res));
    }
    else if (CA_FIELD_IS_NF(CA_FIELD(x, ctx)))
    {
        fmpq_poly_t p;
        fmpq_poly_init(p);

        fmpq_poly_randtest(p, state, qqbar_degree(CA_FIELD_NF_QQBAR(CA_FIELD(x, ctx))), bits);
        fmpz_randtest_not_zero(fmpq_poly_denref(p), state, den_bits);
        fmpz_abs(fmpq_poly_denref(p), fmpq_poly_denref(p));
        fmpq_poly_canonicalise(p);

        ca_set(res, x, ctx);
        nf_elem_set_fmpq_poly(CA_NF_ELEM(res), p, CA_FIELD_NF(CA_FIELD(x, ctx)));
        ca_condense_field(res, ctx);

        fmpq_poly_clear(p);
    }
    else
    {
        flint_throw(FLINT_ERROR, "ca_randtest_same_nf: not implemented\n");
    }
}

void
ca_randtest_rational(ca_t res, flint_rand_t state, slong bits, ca_ctx_t ctx)
{
    fmpq_t t;
    fmpq_init(t);
    fmpq_randtest(t, state, bits);
    ca_set_fmpq(res, t, ctx);
    fmpq_clear(t);
}

void
ca_randtest_qqbar(ca_t res, flint_rand_t state, slong deg, slong bits, ca_ctx_t ctx)
{
    qqbar_t q;
    qqbar_init(q);
    qqbar_randtest(q, state, deg, bits);
    ca_set_qqbar(res, q, ctx);
    qqbar_clear(q);
}

void
ca_randtest(ca_t res, flint_rand_t state, slong len, slong bits, ca_ctx_t ctx)
{
    if (len == 0 || n_randint(state, 2))
    {
        switch (n_randint(state, 10))
        {
            case 0:
                ca_randtest_qqbar(res, state, 2, FLINT_MAX(bits, 2), ctx);
                break;
            case 1:
                ca_i(res, ctx);
                break;
            case 2:
                ca_pi(res, ctx);
                break;
            case 3:
                ca_pi_i(res, ctx);
                break;
            case 4:
                ca_set_si(res, -6 + (slong) n_randint(state, 13), ctx);
                ca_sqrt(res, res, ctx);
                break;
            default:
                ca_randtest_rational(res, state, bits, ctx);
                break;
        }
    }
    else if (n_randint(state, 10) != 0)
    {
        ca_t t;
        ca_init(t, ctx);
        ca_randtest(t, state, len - 1, bits, ctx);
        ca_randtest(res, state, len - 1, bits, ctx);

        switch (n_randint(state, 4))
        {
            case 0:
                ca_add(res, res, t, ctx);
                break;
            case 1:
                ca_sub(res, res, t, ctx);
                break;
            case 2:
                ca_mul(res, res, t, ctx);
                break;
            default:
                if (ca_check_is_zero(t, ctx) == T_FALSE)
                    ca_div(res, res, t, ctx);
                break;
        }

        ca_clear(t, ctx);
    }
    else
    {
        ca_randtest(res, state, len - 1, bits, ctx);

        switch (n_randint(state, 2))
        {
            case 0:
                ca_exp(res, res, ctx);
                break;
            default:
                if (ca_check_is_zero(res, ctx) != T_FALSE)
                    ca_set_ui(res, 2, ctx);
                ca_log(res, res, ctx);
                break;
        }
    }
}

void
ca_randtest_special(ca_t res, flint_rand_t state, slong len, slong bits, ca_ctx_t ctx)
{
    if (n_randint(state, 4) == 0)
    {
        switch (n_randint(state, 7))
        {
            case 0:
                ca_pos_inf(res, ctx);
                break;
            case 1:
                ca_neg_inf(res, ctx);
                break;
            case 2:
                ca_pos_i_inf(res, ctx);
                break;
            case 3:
                ca_neg_i_inf(res, ctx);  /* todo: other signed infs? */
                break;
            case 4:
                ca_uinf(res, ctx);
                break;
            case 5:
                ca_undefined(res, ctx);
                break;
            default:
                ca_unknown(res, ctx);
        }
    }
    else
    {
        ca_randtest(res, state, len, bits, ctx);
    }
}
