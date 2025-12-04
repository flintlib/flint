/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca.h"

void
ca_elliptic_k(ca_t res, const ca_t m, ca_ctx_t ctx)
{
    if (CA_IS_SPECIAL(m))
    {
        if (ca_check_is_pos_inf(m, ctx) == T_TRUE || ca_check_is_neg_inf(m, ctx) == T_TRUE ||
            ca_check_is_pos_i_inf(m, ctx) == T_TRUE || ca_check_is_neg_i_inf(m, ctx) == T_TRUE)
            ca_zero(res, ctx);
        else if (ca_check_is_undefined(m, ctx) == T_TRUE || ca_check_is_uinf(m, ctx) == T_TRUE)
            ca_undefined(res, ctx);
        else
            ca_unknown(res, ctx);
        return;
    }

    if (ca_check_is_zero(m, ctx) == T_TRUE)
    {
        ca_pi(res, ctx);
        ca_div_ui(res, res, 2, ctx);
        return;
    }
    else if (ca_check_is_one(m, ctx) == T_TRUE)
    {
        ca_uinf(res, ctx);
        return;
    }

    _ca_make_field_element(res, _ca_ctx_get_field_fx(ctx, CA_EllipticK, m), ctx);
    fmpz_mpoly_q_gen(CA_MPOLY_Q(res), 0, CA_MCTX_1(ctx));
}

void
ca_elliptic_e(ca_t res, const ca_t m, ca_ctx_t ctx)
{
    if (CA_IS_SPECIAL(m))
    {
        if (ca_check_is_pos_inf(m, ctx) == T_TRUE)
            ca_pos_i_inf(res, ctx);
        else if (ca_check_is_neg_inf(m, ctx) == T_TRUE)
            ca_pos_inf(res, ctx);
        else if (ca_check_is_pos_i_inf(m, ctx) == T_TRUE)
        {
            ca_pos_inf(res, ctx);
            ca_t phase;
            ca_init(phase, ctx);
            ca_i(phase, ctx);
            ca_sqrt(phase, phase, ctx);
            ca_pow_ui(phase, phase, 3, ctx);
            ca_neg(phase, phase, ctx);
            ca_mul(res, res, phase, ctx);
            ca_clear(phase, ctx);
        }
        else if (ca_check_is_neg_i_inf(m, ctx) == T_TRUE)
        {
            ca_pos_inf(res, ctx);
            ca_t phase;
            ca_init(phase, ctx);
            ca_i(phase, ctx);
            ca_sqrt(phase, phase, ctx);
            ca_mul(res, res, phase, ctx);
            ca_clear(phase, ctx);
        }
        else if (ca_check_is_undefined(m, ctx) == T_TRUE || ca_check_is_uinf(m, ctx) == T_TRUE)
            ca_undefined(res, ctx);
        else
            ca_unknown(res, ctx);
        return;
    }

    if (ca_check_is_zero(m, ctx) == T_TRUE)
    {
        ca_pi(res, ctx);
        ca_div_ui(res, res, 2, ctx);
        return;
    }
    else if (ca_check_is_one(m, ctx) == T_TRUE)
    {
        ca_one(res, ctx);
        return;
    }

    _ca_make_field_element(res, _ca_ctx_get_field_fx(ctx, CA_EllipticE, m), ctx);
    fmpz_mpoly_q_gen(CA_MPOLY_Q(res), 0, CA_MCTX_1(ctx));
}

void
ca_elliptic_pi(ca_t res, const ca_t n, const ca_t m, ca_ctx_t ctx)
{
    if (CA_IS_SPECIAL(n))
    {
        if (ca_check_is_pos_inf(n, ctx) == T_TRUE || ca_check_is_neg_inf(n, ctx) == T_TRUE ||
            ca_check_is_pos_i_inf(n, ctx) == T_TRUE || ca_check_is_neg_i_inf(n, ctx) == T_TRUE)
            ca_zero(res, ctx);
        else if (ca_check_is_undefined(n, ctx) == T_TRUE || ca_check_is_uinf(n, ctx) == T_TRUE)
            ca_undefined(res, ctx);
        else
            ca_unknown(res, ctx);
        return;
    }
    
    if (CA_IS_SPECIAL(m))
    {
        if (ca_check_is_pos_inf(m, ctx) == T_TRUE || ca_check_is_neg_inf(m, ctx) == T_TRUE ||
            ca_check_is_pos_i_inf(m, ctx) == T_TRUE || ca_check_is_neg_i_inf(m, ctx) == T_TRUE)
            ca_zero(res, ctx);
        else if (ca_check_is_undefined(m, ctx) == T_TRUE || ca_check_is_uinf(m, ctx) == T_TRUE)
            ca_undefined(res, ctx);
        else
            ca_unknown(res, ctx);
        return;
    }
    
    if (ca_check_is_one(n, ctx) == T_TRUE || ca_check_is_one(m, ctx) == T_TRUE)
    {
        ca_uinf(res, ctx);
        return;
    }

    if (ca_check_is_zero(n, ctx) == T_TRUE)
    {
        _ca_make_field_element(res, _ca_ctx_get_field_fx(ctx, CA_EllipticK, m), ctx);
        fmpz_mpoly_q_gen(CA_MPOLY_Q(res), 0, CA_MCTX_1(ctx));
        return;
    }
    else if (ca_check_is_zero(m, ctx) == T_TRUE)
    {
        ca_t nc;
        ca_init(nc, ctx);
        ca_one(nc, ctx);
        ca_sub(nc, nc, n, ctx);
        ca_sqrt(nc, nc, ctx);
        ca_inv(nc, nc, ctx);
        ca_div_ui(nc, nc, 2, ctx);
        ca_pi(res, ctx);
        ca_mul(res, res, nc, ctx);
        ca_clear(nc, ctx);
        return;
    }

    _ca_make_field_element(res, _ca_ctx_get_field_fxy(ctx, CA_EllipticPi, n, m), ctx);
    fmpz_mpoly_q_gen(CA_MPOLY_Q(res), 0, CA_MCTX_1(ctx));
}

