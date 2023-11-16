/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ca.h"

TEST_FUNCTION_START(ca_properties, state)
{
    {
        ca_ctx_t ctx;
        ca_t x, y;

        ca_ctx_init(ctx);
        ca_init(x, ctx);
        ca_init(y, ctx);

        ca_zero(x, ctx);
        CA_TEST_PROPERTY(ca_check_is_number, "is_number", x, ctx, T_TRUE);
        CA_TEST_PROPERTY(ca_check_is_zero, "is_zero", x, ctx, T_TRUE);
        CA_TEST_PROPERTY(ca_check_is_one, "is_one", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_neg_one, "is_neg_one", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_i, "is_i", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_neg_i, "is_neg_i", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_algebraic, "is_algebraic", x, ctx, T_TRUE);
        CA_TEST_PROPERTY(ca_check_is_rational, "is_rational", x, ctx, T_TRUE);
        CA_TEST_PROPERTY(ca_check_is_integer, "is_integer", x, ctx, T_TRUE);
        CA_TEST_PROPERTY(ca_check_is_real, "is_real", x, ctx, T_TRUE);
        CA_TEST_PROPERTY(ca_check_is_negative_real, "is_negative_real", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_imaginary, "is_imaginary", x, ctx, T_TRUE);
        CA_TEST_PROPERTY(ca_check_is_undefined, "is_undefined", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_infinity, "is_infinity", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_uinf, "is_uinf", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_signed_inf, "is_signed_inf", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_pos_inf, "is_pos_inf", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_neg_inf, "is_neg_inf", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_pos_i_inf, "is_pos_i_inf", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_neg_i_inf, "is_neg_i_inf", x, ctx, T_FALSE);

        ca_one(x, ctx);
        CA_TEST_PROPERTY(ca_check_is_number, "is_number", x, ctx, T_TRUE);
        CA_TEST_PROPERTY(ca_check_is_zero, "is_zero", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_one, "is_one", x, ctx, T_TRUE);
        CA_TEST_PROPERTY(ca_check_is_neg_one, "is_neg_one", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_i, "is_i", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_neg_i, "is_neg_i", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_algebraic, "is_algebraic", x, ctx, T_TRUE);
        CA_TEST_PROPERTY(ca_check_is_rational, "is_rational", x, ctx, T_TRUE);
        CA_TEST_PROPERTY(ca_check_is_integer, "is_integer", x, ctx, T_TRUE);
        CA_TEST_PROPERTY(ca_check_is_real, "is_real", x, ctx, T_TRUE);
        CA_TEST_PROPERTY(ca_check_is_negative_real, "is_negative_real", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_imaginary, "is_imaginary", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_undefined, "is_undefined", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_infinity, "is_infinity", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_uinf, "is_uinf", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_signed_inf, "is_signed_inf", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_pos_inf, "is_pos_inf", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_neg_inf, "is_neg_inf", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_pos_i_inf, "is_pos_i_inf", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_neg_i_inf, "is_neg_i_inf", x, ctx, T_FALSE);

        ca_set_si(x, -1, ctx);
        CA_TEST_PROPERTY(ca_check_is_number, "is_number", x, ctx, T_TRUE);
        CA_TEST_PROPERTY(ca_check_is_zero, "is_zero", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_one, "is_one", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_neg_one, "is_neg_one", x, ctx, T_TRUE);
        CA_TEST_PROPERTY(ca_check_is_i, "is_i", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_neg_i, "is_neg_i", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_algebraic, "is_algebraic", x, ctx, T_TRUE);
        CA_TEST_PROPERTY(ca_check_is_rational, "is_rational", x, ctx, T_TRUE);
        CA_TEST_PROPERTY(ca_check_is_integer, "is_integer", x, ctx, T_TRUE);
        CA_TEST_PROPERTY(ca_check_is_real, "is_real", x, ctx, T_TRUE);
        CA_TEST_PROPERTY(ca_check_is_negative_real, "is_negative_real", x, ctx, T_TRUE);
        CA_TEST_PROPERTY(ca_check_is_imaginary, "is_imaginary", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_undefined, "is_undefined", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_infinity, "is_infinity", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_uinf, "is_uinf", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_signed_inf, "is_signed_inf", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_pos_inf, "is_pos_inf", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_neg_inf, "is_neg_inf", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_pos_i_inf, "is_pos_i_inf", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_neg_i_inf, "is_neg_i_inf", x, ctx, T_FALSE);

        ca_i(x, ctx);
        CA_TEST_PROPERTY(ca_check_is_number, "is_number", x, ctx, T_TRUE);
        CA_TEST_PROPERTY(ca_check_is_zero, "is_zero", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_one, "is_one", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_neg_one, "is_neg_one", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_i, "is_i", x, ctx, T_TRUE);
        CA_TEST_PROPERTY(ca_check_is_neg_i, "is_neg_i", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_algebraic, "is_algebraic", x, ctx, T_TRUE);
        CA_TEST_PROPERTY(ca_check_is_rational, "is_rational", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_integer, "is_integer", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_real, "is_real", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_negative_real, "is_negative_real", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_imaginary, "is_imaginary", x, ctx, T_TRUE);
        CA_TEST_PROPERTY(ca_check_is_undefined, "is_undefined", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_infinity, "is_infinity", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_uinf, "is_uinf", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_signed_inf, "is_signed_inf", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_pos_inf, "is_pos_inf", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_neg_inf, "is_neg_inf", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_pos_i_inf, "is_pos_i_inf", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_neg_i_inf, "is_neg_i_inf", x, ctx, T_FALSE);

        ca_neg_i(x, ctx);
        CA_TEST_PROPERTY(ca_check_is_number, "is_number", x, ctx, T_TRUE);
        CA_TEST_PROPERTY(ca_check_is_zero, "is_zero", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_one, "is_one", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_neg_one, "is_neg_one", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_i, "is_i", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_neg_i, "is_neg_i", x, ctx, T_TRUE);
        CA_TEST_PROPERTY(ca_check_is_algebraic, "is_algebraic", x, ctx, T_TRUE);
        CA_TEST_PROPERTY(ca_check_is_rational, "is_rational", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_integer, "is_integer", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_real, "is_real", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_negative_real, "is_negative_real", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_imaginary, "is_imaginary", x, ctx, T_TRUE);
        CA_TEST_PROPERTY(ca_check_is_undefined, "is_undefined", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_infinity, "is_infinity", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_uinf, "is_uinf", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_signed_inf, "is_signed_inf", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_pos_inf, "is_pos_inf", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_neg_inf, "is_neg_inf", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_pos_i_inf, "is_pos_i_inf", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_neg_i_inf, "is_neg_i_inf", x, ctx, T_FALSE);

        {
            fmpq_t q;
            fmpq_init(q);
            fmpq_set_si(q, -2, 3);
            ca_set_fmpq(x, q, ctx);
            fmpq_clear(q);
        }
        CA_TEST_PROPERTY(ca_check_is_number, "is_number", x, ctx, T_TRUE);
        CA_TEST_PROPERTY(ca_check_is_zero, "is_zero", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_one, "is_one", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_neg_one, "is_neg_one", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_i, "is_i", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_neg_i, "is_neg_i", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_algebraic, "is_algebraic", x, ctx, T_TRUE);
        CA_TEST_PROPERTY(ca_check_is_rational, "is_rational", x, ctx, T_TRUE);
        CA_TEST_PROPERTY(ca_check_is_integer, "is_integer", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_real, "is_real", x, ctx, T_TRUE);
        CA_TEST_PROPERTY(ca_check_is_negative_real, "is_negative_real", x, ctx, T_TRUE);
        CA_TEST_PROPERTY(ca_check_is_imaginary, "is_imaginary", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_undefined, "is_undefined", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_infinity, "is_infinity", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_uinf, "is_uinf", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_signed_inf, "is_signed_inf", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_pos_inf, "is_pos_inf", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_neg_inf, "is_neg_inf", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_pos_i_inf, "is_pos_i_inf", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_neg_i_inf, "is_neg_i_inf", x, ctx, T_FALSE);

        ca_unknown(x, ctx);
        CA_TEST_PROPERTY(ca_check_is_number, "is_number", x, ctx, T_UNKNOWN);
        CA_TEST_PROPERTY(ca_check_is_zero, "is_zero", x, ctx, T_UNKNOWN);
        CA_TEST_PROPERTY(ca_check_is_one, "is_one", x, ctx, T_UNKNOWN);
        CA_TEST_PROPERTY(ca_check_is_neg_one, "is_neg_one", x, ctx, T_UNKNOWN);
        CA_TEST_PROPERTY(ca_check_is_i, "is_i", x, ctx, T_UNKNOWN);
        CA_TEST_PROPERTY(ca_check_is_neg_i, "is_neg_i", x, ctx, T_UNKNOWN);
        CA_TEST_PROPERTY(ca_check_is_algebraic, "is_algebraic", x, ctx, T_UNKNOWN);
        CA_TEST_PROPERTY(ca_check_is_rational, "is_rational", x, ctx, T_UNKNOWN);
        CA_TEST_PROPERTY(ca_check_is_integer, "is_integer", x, ctx, T_UNKNOWN);
        CA_TEST_PROPERTY(ca_check_is_real, "is_real", x, ctx, T_UNKNOWN);
        CA_TEST_PROPERTY(ca_check_is_negative_real, "is_negative_real", x, ctx, T_UNKNOWN);
        CA_TEST_PROPERTY(ca_check_is_imaginary, "is_imaginary", x, ctx, T_UNKNOWN);
        CA_TEST_PROPERTY(ca_check_is_undefined, "is_undefined", x, ctx, T_UNKNOWN);
        CA_TEST_PROPERTY(ca_check_is_infinity, "is_infinity", x, ctx, T_UNKNOWN);
        CA_TEST_PROPERTY(ca_check_is_uinf, "is_uinf", x, ctx, T_UNKNOWN);
        CA_TEST_PROPERTY(ca_check_is_signed_inf, "is_signed_inf", x, ctx, T_UNKNOWN);
        CA_TEST_PROPERTY(ca_check_is_pos_inf, "is_pos_inf", x, ctx, T_UNKNOWN);
        CA_TEST_PROPERTY(ca_check_is_neg_inf, "is_neg_inf", x, ctx, T_UNKNOWN);
        CA_TEST_PROPERTY(ca_check_is_pos_i_inf, "is_pos_i_inf", x, ctx, T_UNKNOWN);
        CA_TEST_PROPERTY(ca_check_is_neg_i_inf, "is_neg_i_inf", x, ctx, T_UNKNOWN);

        ca_undefined(x, ctx);
        CA_TEST_PROPERTY(ca_check_is_number, "is_number", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_zero, "is_zero", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_one, "is_one", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_neg_one, "is_neg_one", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_i, "is_i", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_neg_i, "is_neg_i", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_algebraic, "is_algebraic", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_rational, "is_rational", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_integer, "is_integer", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_real, "is_real", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_negative_real, "is_negative_real", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_imaginary, "is_imaginary", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_undefined, "is_undefined", x, ctx, T_TRUE);
        CA_TEST_PROPERTY(ca_check_is_infinity, "is_infinity", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_uinf, "is_uinf", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_signed_inf, "is_signed_inf", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_pos_inf, "is_pos_inf", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_neg_inf, "is_neg_inf", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_pos_i_inf, "is_pos_i_inf", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_neg_i_inf, "is_neg_i_inf", x, ctx, T_FALSE);

        ca_uinf(x, ctx);
        CA_TEST_PROPERTY(ca_check_is_number, "is_number", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_zero, "is_zero", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_one, "is_one", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_neg_one, "is_neg_one", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_i, "is_i", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_neg_i, "is_neg_i", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_algebraic, "is_algebraic", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_rational, "is_rational", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_integer, "is_integer", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_real, "is_real", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_negative_real, "is_negative_real", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_imaginary, "is_imaginary", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_undefined, "is_undefined", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_infinity, "is_infinity", x, ctx, T_TRUE);
        CA_TEST_PROPERTY(ca_check_is_uinf, "is_uinf", x, ctx, T_TRUE);
        CA_TEST_PROPERTY(ca_check_is_signed_inf, "is_signed_inf", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_pos_inf, "is_pos_inf", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_neg_inf, "is_neg_inf", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_pos_i_inf, "is_pos_i_inf", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_neg_i_inf, "is_neg_i_inf", x, ctx, T_FALSE);

        ca_pos_inf(x, ctx);
        CA_TEST_PROPERTY(ca_check_is_number, "is_number", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_zero, "is_zero", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_one, "is_one", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_neg_one, "is_neg_one", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_i, "is_i", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_neg_i, "is_neg_i", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_algebraic, "is_algebraic", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_rational, "is_rational", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_integer, "is_integer", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_real, "is_real", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_negative_real, "is_negative_real", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_imaginary, "is_imaginary", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_undefined, "is_undefined", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_infinity, "is_infinity", x, ctx, T_TRUE);
        CA_TEST_PROPERTY(ca_check_is_uinf, "is_uinf", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_signed_inf, "is_signed_inf", x, ctx, T_TRUE);
        CA_TEST_PROPERTY(ca_check_is_pos_inf, "is_pos_inf", x, ctx, T_TRUE);
        CA_TEST_PROPERTY(ca_check_is_neg_inf, "is_neg_inf", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_pos_i_inf, "is_pos_i_inf", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_neg_i_inf, "is_neg_i_inf", x, ctx, T_FALSE);

        ca_neg_inf(x, ctx);
        CA_TEST_PROPERTY(ca_check_is_number, "is_number", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_zero, "is_zero", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_one, "is_one", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_neg_one, "is_neg_one", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_i, "is_i", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_neg_i, "is_neg_i", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_algebraic, "is_algebraic", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_rational, "is_rational", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_integer, "is_integer", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_real, "is_real", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_negative_real, "is_negative_real", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_imaginary, "is_imaginary", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_undefined, "is_undefined", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_infinity, "is_infinity", x, ctx, T_TRUE);
        CA_TEST_PROPERTY(ca_check_is_uinf, "is_uinf", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_signed_inf, "is_signed_inf", x, ctx, T_TRUE);
        CA_TEST_PROPERTY(ca_check_is_pos_inf, "is_pos_inf", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_neg_inf, "is_neg_inf", x, ctx, T_TRUE);
        CA_TEST_PROPERTY(ca_check_is_pos_i_inf, "is_pos_i_inf", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_neg_i_inf, "is_neg_i_inf", x, ctx, T_FALSE);

        ca_pos_i_inf(x, ctx);
        CA_TEST_PROPERTY(ca_check_is_number, "is_number", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_zero, "is_zero", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_one, "is_one", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_neg_one, "is_neg_one", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_i, "is_i", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_neg_i, "is_neg_i", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_algebraic, "is_algebraic", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_rational, "is_rational", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_integer, "is_integer", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_real, "is_real", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_negative_real, "is_negative_real", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_imaginary, "is_imaginary", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_undefined, "is_undefined", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_infinity, "is_infinity", x, ctx, T_TRUE);
        CA_TEST_PROPERTY(ca_check_is_uinf, "is_uinf", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_signed_inf, "is_signed_inf", x, ctx, T_TRUE);
        CA_TEST_PROPERTY(ca_check_is_pos_inf, "is_pos_inf", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_neg_inf, "is_neg_inf", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_pos_i_inf, "is_pos_i_inf", x, ctx, T_TRUE);
        CA_TEST_PROPERTY(ca_check_is_neg_i_inf, "is_neg_i_inf", x, ctx, T_FALSE);

        ca_neg_i_inf(x, ctx);
        CA_TEST_PROPERTY(ca_check_is_number, "is_number", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_zero, "is_zero", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_one, "is_one", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_neg_one, "is_neg_one", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_i, "is_i", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_neg_i, "is_neg_i", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_algebraic, "is_algebraic", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_rational, "is_rational", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_integer, "is_integer", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_real, "is_real", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_negative_real, "is_negative_real", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_imaginary, "is_imaginary", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_undefined, "is_undefined", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_infinity, "is_infinity", x, ctx, T_TRUE);
        CA_TEST_PROPERTY(ca_check_is_uinf, "is_uinf", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_signed_inf, "is_signed_inf", x, ctx, T_TRUE);
        CA_TEST_PROPERTY(ca_check_is_pos_inf, "is_pos_inf", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_neg_inf, "is_neg_inf", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_pos_i_inf, "is_pos_i_inf", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_neg_i_inf, "is_neg_i_inf", x, ctx, T_TRUE);

        ca_pi(x, ctx);
        CA_TEST_PROPERTY(ca_check_is_number, "is_number", x, ctx, T_TRUE);
        CA_TEST_PROPERTY(ca_check_is_zero, "is_zero", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_one, "is_one", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_neg_one, "is_neg_one", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_i, "is_i", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_neg_i, "is_neg_i", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_algebraic, "is_algebraic", x, ctx, T_UNKNOWN);   /* todo */
        CA_TEST_PROPERTY(ca_check_is_rational, "is_rational", x, ctx, T_UNKNOWN);  /* todo */
        CA_TEST_PROPERTY(ca_check_is_integer, "is_integer", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_real, "is_real", x, ctx, T_TRUE);
        CA_TEST_PROPERTY(ca_check_is_negative_real, "is_negative_real", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_imaginary, "is_imaginary", x, ctx, T_FALSE);

        ca_pi_i(x, ctx);
        CA_TEST_PROPERTY(ca_check_is_number, "is_number", x, ctx, T_TRUE);
        CA_TEST_PROPERTY(ca_check_is_zero, "is_zero", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_one, "is_one", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_neg_one, "is_neg_one", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_i, "is_i", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_neg_i, "is_neg_i", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_algebraic, "is_algebraic", x, ctx, T_UNKNOWN);   /* todo */
        CA_TEST_PROPERTY(ca_check_is_rational, "is_rational", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_integer, "is_integer", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_real, "is_real", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_negative_real, "is_negative_real", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_imaginary, "is_imaginary", x, ctx, T_TRUE);

        ca_set_si(x, -400, ctx);
        ca_exp(x, x, ctx);
        ca_exp(x, x, ctx);
        CA_TEST_PROPERTY(ca_check_is_number, "is_number", x, ctx, T_TRUE);
        CA_TEST_PROPERTY(ca_check_is_zero, "is_zero", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_one, "is_one", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_neg_one, "is_neg_one", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_i, "is_i", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_neg_i, "is_neg_i", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_algebraic, "is_algebraic", x, ctx, T_UNKNOWN);
        CA_TEST_PROPERTY(ca_check_is_rational, "is_rational", x, ctx, T_UNKNOWN);
        CA_TEST_PROPERTY(ca_check_is_integer, "is_integer", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_real, "is_real", x, ctx, T_TRUE);
        CA_TEST_PROPERTY(ca_check_is_negative_real, "is_negative_real", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_imaginary, "is_imaginary", x, ctx, T_FALSE);

        ca_set_si(x, -1000000, ctx);
        ca_exp(x, x, ctx);
        ca_exp(x, x, ctx);
        CA_TEST_PROPERTY(ca_check_is_number, "is_number", x, ctx, T_TRUE);
        CA_TEST_PROPERTY(ca_check_is_zero, "is_zero", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_one, "is_one", x, ctx, T_UNKNOWN);
        CA_TEST_PROPERTY(ca_check_is_neg_one, "is_neg_one", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_i, "is_i", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_neg_i, "is_neg_i", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_algebraic, "is_algebraic", x, ctx, T_UNKNOWN);
        CA_TEST_PROPERTY(ca_check_is_rational, "is_rational", x, ctx, T_UNKNOWN);
        CA_TEST_PROPERTY(ca_check_is_integer, "is_integer", x, ctx, T_UNKNOWN);
        CA_TEST_PROPERTY(ca_check_is_real, "is_real", x, ctx, T_TRUE);
        CA_TEST_PROPERTY(ca_check_is_negative_real, "is_negative_real", x, ctx, T_FALSE);
        CA_TEST_PROPERTY(ca_check_is_imaginary, "is_imaginary", x, ctx, T_FALSE);

        ca_clear(x, ctx);
        ca_clear(y, ctx);
        ca_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
