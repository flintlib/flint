/* This file is public domain. Author: Fredrik Johansson. */

#include "flint/profiler.h"
#include "ca.h"
#include "ca_ext.h"
#include "ca_field.h"

#define START(expr) flint_printf(">>> "); flint_printf(expr)
#define OUT \
    flint_printf("\n"); \
    ca_print(x, ctx); flint_printf("\n\n");
#define OUT2 \
    flint_printf("\n"); \
    ca_print(x, ctx); flint_printf("\n"); \
    flint_printf(">>> Is zero?\n"); truth_print(ca_check_is_zero(x, ctx)); \
    flint_printf("\n\n");

void
ca_set_qqi_si(ca_t res, slong a, slong b, slong c, slong d, ca_ctx_t ctx)
{
    ca_t t;
    ca_init(t, ctx);
    ca_set_si(t, a, ctx);
    ca_div_si(t, t, b, ctx);
    ca_i(res, ctx);
    ca_mul_si(res, res, c, ctx);
    ca_div_si(res, res, d, ctx);
    ca_add(res, res, t, ctx);
    ca_clear(t, ctx);
}

int main(int argc, char *argv[])
{
    ca_ctx_t ctx;
    ca_t x, y, z;

    TIMEIT_ONCE_START
    ca_ctx_init(ctx);

    ca_init(x, ctx);
    ca_init(y, ctx);
    ca_init(z, ctx);

    START("Exp(Pi*I) + 1");
    ca_pi_i(x, ctx);
    ca_exp(x, x, ctx);
    ca_add_ui(x, x, 1, ctx);
    OUT

    START("Log(-1) / (Pi*I)");
    ca_set_si(x, -1, ctx);
    ca_log(x, x, ctx);
    ca_pi_i(y, ctx);
    ca_div(x, x, y, ctx);
    OUT

    START("Log(-I) / (Pi*I)");
    ca_neg_i(x, ctx);
    ca_log(x, x, ctx);
    ca_pi_i(y, ctx);
    ca_div(x, x, y, ctx);
    OUT

    START("Log(1 / 10^123) / Log(100)");
    ca_set_ui(x, 10, ctx);
    ca_pow_ui(x, x, 123, ctx);
    ca_inv(x, x, ctx);
    ca_log(x, x, ctx);
    ca_set_ui(y, 100, ctx);
    ca_log(y, y, ctx);
    ca_div(x, x, y, ctx);
    OUT

    START("Log(1 + Sqrt(2)) / Log(3 + 2*Sqrt(2))");
    ca_sqrt_ui(x, 2, ctx);
    ca_add_ui(x, x, 1, ctx);
    ca_log(x, x, ctx);
    ca_sqrt_ui(y, 2, ctx);
    ca_mul_ui(y, y, 2, ctx);
    ca_add_ui(y, y, 3, ctx);
    ca_log(y, y, ctx);
    ca_div(x, x, y, ctx);
    OUT

    START("Sqrt(2)*Sqrt(3) - Sqrt(6)");
    ca_sqrt_ui(x, 2, ctx);
    ca_sqrt_ui(y, 3, ctx);
    ca_sqrt_ui(z, 6, ctx);
    ca_mul(x, x, y, ctx);
    ca_sub(x, x, z, ctx);
    OUT

    START("Exp(1+Sqrt(2)) * Exp(1-Sqrt(2)) / (Exp(1)^2)");
    ca_sqrt_ui(x, 2, ctx);
    ca_add_ui(x, x, 1, ctx);
    ca_exp(x, x, ctx);
    ca_sqrt_ui(y, 2, ctx);
    ca_ui_sub(y, 1, y, ctx);
    ca_exp(y, y, ctx);
    ca_one(z, ctx);
    ca_exp(z, z, ctx);
    ca_pow_ui(z, z, 2, ctx);
    ca_mul(x, x, y, ctx);
    ca_div(x, x, z, ctx);
    OUT

    START("I^I - Exp(-Pi/2)");
    ca_i(x, ctx);
    ca_pow(x, x, x, ctx);
    ca_pi(y, ctx);
    ca_div_ui(y, y, 2, ctx);
    ca_neg(y, y, ctx);
    ca_exp(y, y, ctx);
    ca_sub(x, x, y, ctx);
    OUT

    START("Exp(Sqrt(3))^2 - Exp(Sqrt(12))");
    ca_sqrt_ui(x, 3, ctx);
    ca_exp(x, x, ctx);
    ca_pow_ui(x, x, 2, ctx);
    ca_sqrt_ui(y, 12, ctx);
    ca_exp(y, y, ctx);
    ca_sub(x, x, y, ctx);
    OUT

    START("2*Log(Pi*I) - 4*Log(Sqrt(Pi)) - Pi*I");
    ca_pi_i(x, ctx);
    ca_log(x, x, ctx);
    ca_mul_ui(x, x, 2, ctx);
    ca_pi(y, ctx);
    ca_sqrt(y, y, ctx);
    ca_log(y, y, ctx);
    ca_mul_ui(y, y, 4, ctx);
    ca_sub(x, x, y, ctx);
    ca_pi_i(y, ctx);
    ca_sub(x, x, y, ctx);
    OUT

    /* Example 1 in [BBK2014] */
    START("-I*Pi/8*Log(2/3-2*I/3)^2 + I*Pi/8*Log(2/3+2*I/3)^2 + Pi^2/12*Log(-1-I) + Pi^2/12*Log(-1+I) + Pi^2/12*Log(1/3-I/3) + Pi^2/12*Log(1/3+I/3) - Pi^2/48*Log(18)");
    ca_zero(x, ctx);
    ca_set_qqi_si(y, 2, 3, -2, 3, ctx);
    ca_log(y, y, ctx);
    ca_sqr(y, y, ctx);
    ca_pi_i(z, ctx);
    ca_mul(y, y, z, ctx);
    ca_div_si(y, y, -8, ctx);
    ca_add(x, x, y, ctx);

    ca_set_qqi_si(y, 2, 3, 2, 3, ctx);
    ca_log(y, y, ctx);
    ca_sqr(y, y, ctx);
    ca_pi_i(z, ctx);
    ca_mul(y, y, z, ctx);
    ca_div_si(y, y, 8, ctx);
    ca_add(x, x, y, ctx);

    ca_set_qqi_si(y, -1, 1, -1, 1, ctx);
    ca_log(y, y, ctx);
    ca_pi(z, ctx);
    ca_sqr(z, z, ctx);
    ca_mul(y, y, z, ctx);
    ca_div_si(y, y, 12, ctx);
    ca_add(x, x, y, ctx);

    ca_set_qqi_si(y, -1, 1, 1, 1, ctx);
    ca_log(y, y, ctx);
    ca_pi(z, ctx);
    ca_sqr(z, z, ctx);
    ca_mul(y, y, z, ctx);
    ca_div_si(y, y, 12, ctx);
    ca_add(x, x, y, ctx);

    ca_set_qqi_si(y, 1, 3, -1, 3, ctx);
    ca_log(y, y, ctx);
    ca_pi(z, ctx);
    ca_sqr(z, z, ctx);
    ca_mul(y, y, z, ctx);
    ca_div_si(y, y, 12, ctx);
    ca_add(x, x, y, ctx);

    ca_set_qqi_si(y, 1, 3, 1, 3, ctx);
    ca_log(y, y, ctx);
    ca_pi(z, ctx);
    ca_sqr(z, z, ctx);
    ca_mul(y, y, z, ctx);
    ca_div_si(y, y, 12, ctx);
    ca_add(x, x, y, ctx);

    ca_set_ui(y, 18, ctx);
    ca_log(y, y, ctx);
    ca_pi(z, ctx);
    ca_sqr(z, z, ctx);
    ca_mul(y, y, z, ctx);
    ca_div_si(y, y, -48, ctx);
    ca_sub(x, x, y, ctx);
    OUT

    START("Sqrt(5 + 2*Sqrt(6)) - Sqrt(2) - Sqrt(3)");
    ca_sqrt_ui(x, 6, ctx);
    ca_mul_ui(x, x, 2, ctx);
    ca_add_ui(x, x, 5, ctx);
    ca_sqrt(x, x, ctx);
    ca_sqrt_ui(y, 2, ctx);
    ca_sub(x, x, y, ctx);
    ca_sqrt_ui(y, 3, ctx);
    ca_sub(x, x, y, ctx);
    OUT2

    START("Sqrt(I) - (1+I)/Sqrt(2)");
    ca_i(x, ctx);
    ca_sqrt(x, x, ctx);
    ca_i(y, ctx);
    ca_add_ui(y, y, 1, ctx);
    ca_sqrt_ui(z, 2, ctx);
    ca_div(y, y, z, ctx);
    ca_sub(x, x, y, ctx);
    OUT2

    START("Exp(Pi*Sqrt(163)) - (640320^3 + 744)");
    ca_pi(x, ctx);
    ca_sqrt_ui(y, 163, ctx);
    ca_mul(x, x, y, ctx);
    ca_exp(x, x, ctx);
    ca_set_ui(y, 640320, ctx);
    ca_pow_ui(y, y, 3, ctx);
    ca_add_ui(y, y, 744, ctx);
    ca_sub(x, x, y, ctx);
    OUT

    /* Taken (slightly tweaked) from: https://reference.wolfram.com/language/ref/PossibleZeroQ.html */
    START("Erf(2*Log(Sqrt(1/2-Sqrt(2)/4))+Log(4)) - Erf(Log(2-Sqrt(2)))");
    ca_sqrt_ui(x, 2, ctx);
    ca_div_ui(x, x, 4, ctx);
    ca_one(y, ctx);
    ca_div_ui(y, y, 2, ctx);
    ca_sub(x, y, x, ctx);
    ca_sqrt(x, x, ctx);
    ca_log(x, x, ctx);
    ca_mul_ui(x, x, 2, ctx);
    ca_set_ui(y, 4, ctx);
    ca_log(y, y, ctx);
    ca_add(x, y, x, ctx);
    ca_erf(x, x, ctx);
    ca_sqrt_ui(y, 2, ctx);
    ca_ui_sub(y, 2, y, ctx);
    ca_log(y, y, ctx);
    ca_erf(y, y, ctx);
    ca_sub(x, x, y, ctx);
    OUT

    flint_printf("\n");
    ca_clear(x, ctx);
    ca_clear(y, ctx);
    ca_clear(z, ctx);
    ca_ctx_clear(ctx);

    TIMEIT_ONCE_STOP
    SHOW_MEMORY_USAGE

    flint_cleanup();
    return EXIT_SUCCESS;
}
