/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca.h"

#define TEST(f, s, x, ctx, expected) \
     do { \
        truth_t t; \
        t = f(x, ctx); \
        if (t != expected) \
        { \
            flint_printf("FAIL\n"); \
            flint_printf("%s\n\n", s); \
            flint_printf("x = "); ca_print(x, ctx); flint_printf("\n\n"); \
            flint_printf("got = "); truth_print(t); flint_printf("\n\n"); \
            flint_printf("expected = "); truth_print(expected); flint_printf("\n\n"); \
            flint_abort(); \
        } \
     } while (0) \

int main()
{
    flint_rand_t state;

    flint_printf("properties....");
    fflush(stdout);

    flint_randinit(state);

    {
        ca_ctx_t ctx;
        ca_t x;

        ca_ctx_init(ctx);
        ca_init(x, ctx);

        ca_zero(x, ctx);
        TEST(ca_check_is_number, "is_number", x, ctx, T_TRUE);
        TEST(ca_check_is_zero, "is_zero", x, ctx, T_TRUE);
        TEST(ca_check_is_one, "is_one", x, ctx, T_FALSE);
        TEST(ca_check_is_neg_one, "is_neg_one", x, ctx, T_FALSE);
        TEST(ca_check_is_i, "is_i", x, ctx, T_FALSE);
        TEST(ca_check_is_neg_i, "is_neg_i", x, ctx, T_FALSE);
        TEST(ca_check_is_algebraic, "is_algebraic", x, ctx, T_TRUE);
        TEST(ca_check_is_rational, "is_rational", x, ctx, T_TRUE);
        TEST(ca_check_is_integer, "is_integer", x, ctx, T_TRUE);
        TEST(ca_check_is_real, "is_real", x, ctx, T_TRUE);
        TEST(ca_check_is_imaginary, "is_imaginary", x, ctx, T_TRUE);
        TEST(ca_check_is_nonreal, "is_nonreal", x, ctx, T_FALSE);
        TEST(ca_check_is_undefined, "is_undefined", x, ctx, T_FALSE);
        TEST(ca_check_is_infinity, "is_infinity", x, ctx, T_FALSE);
        TEST(ca_check_is_uinf, "is_uinf", x, ctx, T_FALSE);
        TEST(ca_check_is_signed_inf, "is_signed_inf", x, ctx, T_FALSE);
        TEST(ca_check_is_pos_inf, "is_pos_inf", x, ctx, T_FALSE);
        TEST(ca_check_is_neg_inf, "is_neg_inf", x, ctx, T_FALSE);
        TEST(ca_check_is_pos_i_inf, "is_pos_i_inf", x, ctx, T_FALSE);
        TEST(ca_check_is_neg_i_inf, "is_neg_i_inf", x, ctx, T_FALSE);

        ca_one(x, ctx);
        TEST(ca_check_is_number, "is_number", x, ctx, T_TRUE);
        TEST(ca_check_is_zero, "is_zero", x, ctx, T_FALSE);
        TEST(ca_check_is_one, "is_one", x, ctx, T_TRUE);
        TEST(ca_check_is_neg_one, "is_neg_one", x, ctx, T_FALSE);
        TEST(ca_check_is_i, "is_i", x, ctx, T_FALSE);
        TEST(ca_check_is_neg_i, "is_neg_i", x, ctx, T_FALSE);
        TEST(ca_check_is_algebraic, "is_algebraic", x, ctx, T_TRUE);
        TEST(ca_check_is_rational, "is_rational", x, ctx, T_TRUE);
        TEST(ca_check_is_integer, "is_integer", x, ctx, T_TRUE);
        TEST(ca_check_is_real, "is_real", x, ctx, T_TRUE);
        TEST(ca_check_is_imaginary, "is_imaginary", x, ctx, T_FALSE);
        TEST(ca_check_is_nonreal, "is_nonreal", x, ctx, T_FALSE);
        TEST(ca_check_is_undefined, "is_undefined", x, ctx, T_FALSE);
        TEST(ca_check_is_infinity, "is_infinity", x, ctx, T_FALSE);
        TEST(ca_check_is_uinf, "is_uinf", x, ctx, T_FALSE);
        TEST(ca_check_is_signed_inf, "is_signed_inf", x, ctx, T_FALSE);
        TEST(ca_check_is_pos_inf, "is_pos_inf", x, ctx, T_FALSE);
        TEST(ca_check_is_neg_inf, "is_neg_inf", x, ctx, T_FALSE);
        TEST(ca_check_is_pos_i_inf, "is_pos_i_inf", x, ctx, T_FALSE);
        TEST(ca_check_is_neg_i_inf, "is_neg_i_inf", x, ctx, T_FALSE);

        ca_set_si(x, -1, ctx);
        TEST(ca_check_is_number, "is_number", x, ctx, T_TRUE);
        TEST(ca_check_is_zero, "is_zero", x, ctx, T_FALSE);
        TEST(ca_check_is_one, "is_one", x, ctx, T_FALSE);
        TEST(ca_check_is_neg_one, "is_neg_one", x, ctx, T_TRUE);
        TEST(ca_check_is_i, "is_i", x, ctx, T_FALSE);
        TEST(ca_check_is_neg_i, "is_neg_i", x, ctx, T_FALSE);
        TEST(ca_check_is_algebraic, "is_algebraic", x, ctx, T_TRUE);
        TEST(ca_check_is_rational, "is_rational", x, ctx, T_TRUE);
        TEST(ca_check_is_integer, "is_integer", x, ctx, T_TRUE);
        TEST(ca_check_is_real, "is_real", x, ctx, T_TRUE);
        TEST(ca_check_is_imaginary, "is_imaginary", x, ctx, T_FALSE);
        TEST(ca_check_is_nonreal, "is_nonreal", x, ctx, T_FALSE);
        TEST(ca_check_is_undefined, "is_undefined", x, ctx, T_FALSE);
        TEST(ca_check_is_infinity, "is_infinity", x, ctx, T_FALSE);
        TEST(ca_check_is_uinf, "is_uinf", x, ctx, T_FALSE);
        TEST(ca_check_is_signed_inf, "is_signed_inf", x, ctx, T_FALSE);
        TEST(ca_check_is_pos_inf, "is_pos_inf", x, ctx, T_FALSE);
        TEST(ca_check_is_neg_inf, "is_neg_inf", x, ctx, T_FALSE);
        TEST(ca_check_is_pos_i_inf, "is_pos_i_inf", x, ctx, T_FALSE);
        TEST(ca_check_is_neg_i_inf, "is_neg_i_inf", x, ctx, T_FALSE);

        ca_i(x, ctx);
        TEST(ca_check_is_number, "is_number", x, ctx, T_TRUE);
        TEST(ca_check_is_zero, "is_zero", x, ctx, T_FALSE);
        TEST(ca_check_is_one, "is_one", x, ctx, T_FALSE);
        TEST(ca_check_is_neg_one, "is_neg_one", x, ctx, T_FALSE);
        TEST(ca_check_is_i, "is_i", x, ctx, T_TRUE);
        TEST(ca_check_is_neg_i, "is_neg_i", x, ctx, T_FALSE);
        TEST(ca_check_is_algebraic, "is_algebraic", x, ctx, T_TRUE);
        TEST(ca_check_is_rational, "is_rational", x, ctx, T_FALSE);
        TEST(ca_check_is_integer, "is_integer", x, ctx, T_FALSE);
        TEST(ca_check_is_real, "is_real", x, ctx, T_FALSE);
        TEST(ca_check_is_imaginary, "is_imaginary", x, ctx, T_TRUE);
        TEST(ca_check_is_nonreal, "is_nonreal", x, ctx, T_TRUE);
        TEST(ca_check_is_undefined, "is_undefined", x, ctx, T_FALSE);
        TEST(ca_check_is_infinity, "is_infinity", x, ctx, T_FALSE);
        TEST(ca_check_is_uinf, "is_uinf", x, ctx, T_FALSE);
        TEST(ca_check_is_signed_inf, "is_signed_inf", x, ctx, T_FALSE);
        TEST(ca_check_is_pos_inf, "is_pos_inf", x, ctx, T_FALSE);
        TEST(ca_check_is_neg_inf, "is_neg_inf", x, ctx, T_FALSE);
        TEST(ca_check_is_pos_i_inf, "is_pos_i_inf", x, ctx, T_FALSE);
        TEST(ca_check_is_neg_i_inf, "is_neg_i_inf", x, ctx, T_FALSE);

/*
        ca_neg_i(x, ctx);
        TEST(ca_check_is_number, "is_number", x, ctx, T_TRUE);
        TEST(ca_check_is_zero, "is_zero", x, ctx, T_FALSE);
        TEST(ca_check_is_one, "is_one", x, ctx, T_FALSE);
        TEST(ca_check_is_neg_one, "is_neg_one", x, ctx, T_FALSE);
        TEST(ca_check_is_i, "is_i", x, ctx, T_TRUE);
        TEST(ca_check_is_neg_i, "is_neg_i", x, ctx, T_FALSE);
        TEST(ca_check_is_algebraic, "is_algebraic", x, ctx, T_TRUE);
        TEST(ca_check_is_rational, "is_rational", x, ctx, T_FALSE);
        TEST(ca_check_is_integer, "is_integer", x, ctx, T_FALSE);
        TEST(ca_check_is_real, "is_real", x, ctx, T_FALSE);
        TEST(ca_check_is_imaginary, "is_imaginary", x, ctx, T_TRUE);
        TEST(ca_check_is_nonreal, "is_nonreal", x, ctx, T_TRUE);
        TEST(ca_check_is_undefined, "is_undefined", x, ctx, T_FALSE);
        TEST(ca_check_is_infinity, "is_infinity", x, ctx, T_FALSE);
        TEST(ca_check_is_uinf, "is_uinf", x, ctx, T_FALSE);
        TEST(ca_check_is_signed_inf, "is_signed_inf", x, ctx, T_FALSE);
        TEST(ca_check_is_pos_inf, "is_pos_inf", x, ctx, T_FALSE);
        TEST(ca_check_is_neg_inf, "is_neg_inf", x, ctx, T_FALSE);
        TEST(ca_check_is_pos_i_inf, "is_pos_i_inf", x, ctx, T_FALSE);
        TEST(ca_check_is_neg_i_inf, "is_neg_i_inf", x, ctx, T_FALSE);
*/

        {
            fmpq_t q;
            fmpq_init(q);
            fmpq_set_si(q, -2, 3);
            ca_set_fmpq(x, q, ctx);
            fmpq_clear(q);
        }
        TEST(ca_check_is_number, "is_number", x, ctx, T_TRUE);
        TEST(ca_check_is_zero, "is_zero", x, ctx, T_FALSE);
        TEST(ca_check_is_one, "is_one", x, ctx, T_FALSE);
        TEST(ca_check_is_neg_one, "is_neg_one", x, ctx, T_FALSE);
        TEST(ca_check_is_i, "is_i", x, ctx, T_FALSE);
        TEST(ca_check_is_neg_i, "is_neg_i", x, ctx, T_FALSE);
        TEST(ca_check_is_algebraic, "is_algebraic", x, ctx, T_TRUE);
        TEST(ca_check_is_rational, "is_rational", x, ctx, T_TRUE);
        TEST(ca_check_is_integer, "is_integer", x, ctx, T_FALSE);
        TEST(ca_check_is_real, "is_real", x, ctx, T_TRUE);
        TEST(ca_check_is_imaginary, "is_imaginary", x, ctx, T_FALSE);
        TEST(ca_check_is_nonreal, "is_nonreal", x, ctx, T_FALSE);
        TEST(ca_check_is_undefined, "is_undefined", x, ctx, T_FALSE);
        TEST(ca_check_is_infinity, "is_infinity", x, ctx, T_FALSE);
        TEST(ca_check_is_uinf, "is_uinf", x, ctx, T_FALSE);
        TEST(ca_check_is_signed_inf, "is_signed_inf", x, ctx, T_FALSE);
        TEST(ca_check_is_pos_inf, "is_pos_inf", x, ctx, T_FALSE);
        TEST(ca_check_is_neg_inf, "is_neg_inf", x, ctx, T_FALSE);
        TEST(ca_check_is_pos_i_inf, "is_pos_i_inf", x, ctx, T_FALSE);
        TEST(ca_check_is_neg_i_inf, "is_neg_i_inf", x, ctx, T_FALSE);

        ca_unknown(x, ctx);
        TEST(ca_check_is_number, "is_number", x, ctx, T_UNKNOWN);
        TEST(ca_check_is_zero, "is_zero", x, ctx, T_UNKNOWN);
        TEST(ca_check_is_one, "is_one", x, ctx, T_UNKNOWN);
        TEST(ca_check_is_neg_one, "is_neg_one", x, ctx, T_UNKNOWN);
        TEST(ca_check_is_i, "is_i", x, ctx, T_UNKNOWN);
        TEST(ca_check_is_neg_i, "is_neg_i", x, ctx, T_UNKNOWN);
        TEST(ca_check_is_algebraic, "is_algebraic", x, ctx, T_UNKNOWN);
        TEST(ca_check_is_rational, "is_rational", x, ctx, T_UNKNOWN);
        TEST(ca_check_is_integer, "is_integer", x, ctx, T_UNKNOWN);
        TEST(ca_check_is_real, "is_real", x, ctx, T_UNKNOWN);
        TEST(ca_check_is_imaginary, "is_imaginary", x, ctx, T_UNKNOWN);
        TEST(ca_check_is_nonreal, "is_nonreal", x, ctx, T_UNKNOWN);
        TEST(ca_check_is_undefined, "is_undefined", x, ctx, T_UNKNOWN);
        TEST(ca_check_is_infinity, "is_infinity", x, ctx, T_UNKNOWN);
        TEST(ca_check_is_uinf, "is_uinf", x, ctx, T_UNKNOWN);
        TEST(ca_check_is_signed_inf, "is_signed_inf", x, ctx, T_UNKNOWN);
        TEST(ca_check_is_pos_inf, "is_pos_inf", x, ctx, T_UNKNOWN);
        TEST(ca_check_is_neg_inf, "is_neg_inf", x, ctx, T_UNKNOWN);
        TEST(ca_check_is_pos_i_inf, "is_pos_i_inf", x, ctx, T_UNKNOWN);
        TEST(ca_check_is_neg_i_inf, "is_neg_i_inf", x, ctx, T_UNKNOWN);

        ca_undefined(x, ctx);
        TEST(ca_check_is_number, "is_number", x, ctx, T_FALSE);
        TEST(ca_check_is_zero, "is_zero", x, ctx, T_FALSE);
        TEST(ca_check_is_one, "is_one", x, ctx, T_FALSE);
        TEST(ca_check_is_neg_one, "is_neg_one", x, ctx, T_FALSE);
        TEST(ca_check_is_i, "is_i", x, ctx, T_FALSE);
        TEST(ca_check_is_neg_i, "is_neg_i", x, ctx, T_FALSE);
        TEST(ca_check_is_algebraic, "is_algebraic", x, ctx, T_FALSE);
        TEST(ca_check_is_rational, "is_rational", x, ctx, T_FALSE);
        TEST(ca_check_is_integer, "is_integer", x, ctx, T_FALSE);
        TEST(ca_check_is_real, "is_real", x, ctx, T_FALSE);
        TEST(ca_check_is_imaginary, "is_imaginary", x, ctx, T_FALSE);
        TEST(ca_check_is_nonreal, "is_nonreal", x, ctx, T_FALSE);
        TEST(ca_check_is_undefined, "is_undefined", x, ctx, T_TRUE);
        TEST(ca_check_is_infinity, "is_infinity", x, ctx, T_FALSE);
        TEST(ca_check_is_uinf, "is_uinf", x, ctx, T_FALSE);
        TEST(ca_check_is_signed_inf, "is_signed_inf", x, ctx, T_FALSE);
        TEST(ca_check_is_pos_inf, "is_pos_inf", x, ctx, T_FALSE);
        TEST(ca_check_is_neg_inf, "is_neg_inf", x, ctx, T_FALSE);
        TEST(ca_check_is_pos_i_inf, "is_pos_i_inf", x, ctx, T_FALSE);
        TEST(ca_check_is_neg_i_inf, "is_neg_i_inf", x, ctx, T_FALSE);

        ca_uinf(x, ctx);
        TEST(ca_check_is_number, "is_number", x, ctx, T_FALSE);
        TEST(ca_check_is_zero, "is_zero", x, ctx, T_FALSE);
        TEST(ca_check_is_one, "is_one", x, ctx, T_FALSE);
        TEST(ca_check_is_neg_one, "is_neg_one", x, ctx, T_FALSE);
        TEST(ca_check_is_i, "is_i", x, ctx, T_FALSE);
        TEST(ca_check_is_neg_i, "is_neg_i", x, ctx, T_FALSE);
        TEST(ca_check_is_algebraic, "is_algebraic", x, ctx, T_FALSE);
        TEST(ca_check_is_rational, "is_rational", x, ctx, T_FALSE);
        TEST(ca_check_is_integer, "is_integer", x, ctx, T_FALSE);
        TEST(ca_check_is_real, "is_real", x, ctx, T_FALSE);
        TEST(ca_check_is_imaginary, "is_imaginary", x, ctx, T_FALSE);
        TEST(ca_check_is_nonreal, "is_nonreal", x, ctx, T_FALSE);
        TEST(ca_check_is_undefined, "is_undefined", x, ctx, T_FALSE);
        TEST(ca_check_is_infinity, "is_infinity", x, ctx, T_TRUE);
        TEST(ca_check_is_uinf, "is_uinf", x, ctx, T_TRUE);
        TEST(ca_check_is_signed_inf, "is_signed_inf", x, ctx, T_FALSE);
        TEST(ca_check_is_pos_inf, "is_pos_inf", x, ctx, T_FALSE);
        TEST(ca_check_is_neg_inf, "is_neg_inf", x, ctx, T_FALSE);
        TEST(ca_check_is_pos_i_inf, "is_pos_i_inf", x, ctx, T_FALSE);
        TEST(ca_check_is_neg_i_inf, "is_neg_i_inf", x, ctx, T_FALSE);

        ca_pos_inf(x, ctx);
        TEST(ca_check_is_number, "is_number", x, ctx, T_FALSE);
        TEST(ca_check_is_zero, "is_zero", x, ctx, T_FALSE);
        TEST(ca_check_is_one, "is_one", x, ctx, T_FALSE);
        TEST(ca_check_is_neg_one, "is_neg_one", x, ctx, T_FALSE);
        TEST(ca_check_is_i, "is_i", x, ctx, T_FALSE);
        TEST(ca_check_is_neg_i, "is_neg_i", x, ctx, T_FALSE);
        TEST(ca_check_is_algebraic, "is_algebraic", x, ctx, T_FALSE);
        TEST(ca_check_is_rational, "is_rational", x, ctx, T_FALSE);
        TEST(ca_check_is_integer, "is_integer", x, ctx, T_FALSE);
        TEST(ca_check_is_real, "is_real", x, ctx, T_FALSE);
        TEST(ca_check_is_imaginary, "is_imaginary", x, ctx, T_FALSE);
        TEST(ca_check_is_nonreal, "is_nonreal", x, ctx, T_FALSE);
        TEST(ca_check_is_undefined, "is_undefined", x, ctx, T_FALSE);
        TEST(ca_check_is_infinity, "is_infinity", x, ctx, T_TRUE);
        TEST(ca_check_is_uinf, "is_uinf", x, ctx, T_FALSE);
        TEST(ca_check_is_signed_inf, "is_signed_inf", x, ctx, T_TRUE);
        TEST(ca_check_is_pos_inf, "is_pos_inf", x, ctx, T_TRUE);
        TEST(ca_check_is_neg_inf, "is_neg_inf", x, ctx, T_FALSE);
        TEST(ca_check_is_pos_i_inf, "is_pos_i_inf", x, ctx, T_FALSE);
        TEST(ca_check_is_neg_i_inf, "is_neg_i_inf", x, ctx, T_FALSE);

        ca_neg_inf(x, ctx);
        TEST(ca_check_is_number, "is_number", x, ctx, T_FALSE);
        TEST(ca_check_is_zero, "is_zero", x, ctx, T_FALSE);
        TEST(ca_check_is_one, "is_one", x, ctx, T_FALSE);
        TEST(ca_check_is_neg_one, "is_neg_one", x, ctx, T_FALSE);
        TEST(ca_check_is_i, "is_i", x, ctx, T_FALSE);
        TEST(ca_check_is_neg_i, "is_neg_i", x, ctx, T_FALSE);
        TEST(ca_check_is_algebraic, "is_algebraic", x, ctx, T_FALSE);
        TEST(ca_check_is_rational, "is_rational", x, ctx, T_FALSE);
        TEST(ca_check_is_integer, "is_integer", x, ctx, T_FALSE);
        TEST(ca_check_is_real, "is_real", x, ctx, T_FALSE);
        TEST(ca_check_is_imaginary, "is_imaginary", x, ctx, T_FALSE);
        TEST(ca_check_is_nonreal, "is_nonreal", x, ctx, T_FALSE);
        TEST(ca_check_is_undefined, "is_undefined", x, ctx, T_FALSE);
        TEST(ca_check_is_infinity, "is_infinity", x, ctx, T_TRUE);
        TEST(ca_check_is_uinf, "is_uinf", x, ctx, T_FALSE);
        TEST(ca_check_is_signed_inf, "is_signed_inf", x, ctx, T_TRUE);
        TEST(ca_check_is_pos_inf, "is_pos_inf", x, ctx, T_FALSE);
        TEST(ca_check_is_neg_inf, "is_neg_inf", x, ctx, T_TRUE);
        TEST(ca_check_is_pos_i_inf, "is_pos_i_inf", x, ctx, T_FALSE);
        TEST(ca_check_is_neg_i_inf, "is_neg_i_inf", x, ctx, T_FALSE);

        ca_pos_i_inf(x, ctx);
        TEST(ca_check_is_number, "is_number", x, ctx, T_FALSE);
        TEST(ca_check_is_zero, "is_zero", x, ctx, T_FALSE);
        TEST(ca_check_is_one, "is_one", x, ctx, T_FALSE);
        TEST(ca_check_is_neg_one, "is_neg_one", x, ctx, T_FALSE);
        TEST(ca_check_is_i, "is_i", x, ctx, T_FALSE);
        TEST(ca_check_is_neg_i, "is_neg_i", x, ctx, T_FALSE);
        TEST(ca_check_is_algebraic, "is_algebraic", x, ctx, T_FALSE);
        TEST(ca_check_is_rational, "is_rational", x, ctx, T_FALSE);
        TEST(ca_check_is_integer, "is_integer", x, ctx, T_FALSE);
        TEST(ca_check_is_real, "is_real", x, ctx, T_FALSE);
        TEST(ca_check_is_imaginary, "is_imaginary", x, ctx, T_FALSE);
        TEST(ca_check_is_nonreal, "is_nonreal", x, ctx, T_FALSE);
        TEST(ca_check_is_undefined, "is_undefined", x, ctx, T_FALSE);
        TEST(ca_check_is_infinity, "is_infinity", x, ctx, T_TRUE);
        TEST(ca_check_is_uinf, "is_uinf", x, ctx, T_FALSE);
        TEST(ca_check_is_signed_inf, "is_signed_inf", x, ctx, T_TRUE);
        TEST(ca_check_is_pos_inf, "is_pos_inf", x, ctx, T_FALSE);
        TEST(ca_check_is_neg_inf, "is_neg_inf", x, ctx, T_FALSE);
        TEST(ca_check_is_pos_i_inf, "is_pos_i_inf", x, ctx, T_TRUE);
        TEST(ca_check_is_neg_i_inf, "is_neg_i_inf", x, ctx, T_FALSE);

        ca_neg_i_inf(x, ctx);
        TEST(ca_check_is_number, "is_number", x, ctx, T_FALSE);
        TEST(ca_check_is_zero, "is_zero", x, ctx, T_FALSE);
        TEST(ca_check_is_one, "is_one", x, ctx, T_FALSE);
        TEST(ca_check_is_neg_one, "is_neg_one", x, ctx, T_FALSE);
        TEST(ca_check_is_i, "is_i", x, ctx, T_FALSE);
        TEST(ca_check_is_neg_i, "is_neg_i", x, ctx, T_FALSE);
        TEST(ca_check_is_algebraic, "is_algebraic", x, ctx, T_FALSE);
        TEST(ca_check_is_rational, "is_rational", x, ctx, T_FALSE);
        TEST(ca_check_is_integer, "is_integer", x, ctx, T_FALSE);
        TEST(ca_check_is_real, "is_real", x, ctx, T_FALSE);
        TEST(ca_check_is_imaginary, "is_imaginary", x, ctx, T_FALSE);
        TEST(ca_check_is_nonreal, "is_nonreal", x, ctx, T_FALSE);
        TEST(ca_check_is_undefined, "is_undefined", x, ctx, T_FALSE);
        TEST(ca_check_is_infinity, "is_infinity", x, ctx, T_TRUE);
        TEST(ca_check_is_uinf, "is_uinf", x, ctx, T_FALSE);
        TEST(ca_check_is_signed_inf, "is_signed_inf", x, ctx, T_TRUE);
        TEST(ca_check_is_pos_inf, "is_pos_inf", x, ctx, T_FALSE);
        TEST(ca_check_is_neg_inf, "is_neg_inf", x, ctx, T_FALSE);
        TEST(ca_check_is_pos_i_inf, "is_pos_i_inf", x, ctx, T_FALSE);
        TEST(ca_check_is_neg_i_inf, "is_neg_i_inf", x, ctx, T_TRUE);

        ca_clear(x, ctx);
        ca_ctx_clear(ctx);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

