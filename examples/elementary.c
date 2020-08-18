/* This file is public domain. Author: Fredrik Johansson. */

#include "flint/profiler.h"
#include "ca.h"
#include "ca_ext.h"
#include "ca_field.h"

int main(int argc, char *argv[])
{
    ca_ctx_t ctx;
    ca_t x, y, z;

    TIMEIT_ONCE_START
    ca_ctx_init(ctx);

    ca_init(x, ctx);
    ca_init(y, ctx);
    ca_init(z, ctx);

    flint_printf("1:  exp(pi*i) + 1                                = ");
    ca_pi_i(x, ctx);
    ca_exp(x, x, ctx);
    ca_add_ui(x, x, 1, ctx);
    ca_print(x, ctx);
    flint_printf("\n");

    flint_printf("2:  log(-1)/(pi*i)                               = ");
    ca_set_si(x, -1, ctx);
    ca_log(x, x, ctx);
    ca_pi_i(y, ctx);
    ca_div(x, x, y, ctx);
    ca_print(x, ctx);
    flint_printf("\n");

    flint_printf("3:  log(-i)/(pi*i)                               = ");
    ca_neg_i(x, ctx);
    ca_log(x, x, ctx);
    ca_pi_i(y, ctx);
    ca_div(x, x, y, ctx);
    ca_print(x, ctx);
    flint_printf("\n");

    flint_printf("4:  log(1/10^123)/log(100)                       = ");
    ca_set_ui(x, 10, ctx);
    ca_pow_ui(x, x, 123, ctx);
    ca_inv(x, x, ctx);
    ca_log(x, x, ctx);
    ca_set_ui(y, 100, ctx);
    ca_log(y, y, ctx);
    ca_div(x, x, y, ctx);
    ca_print(x, ctx);
    flint_printf("\n");

    flint_printf("5:  log(1+sqrt(2)) / log(3+2*sqrt(2))            = ");
    ca_sqrt_ui(x, 2, ctx);
    ca_add_ui(x, x, 1, ctx);
    ca_log(x, x, ctx);
    ca_sqrt_ui(y, 2, ctx);
    ca_mul_ui(y, y, 2, ctx);
    ca_add_ui(y, y, 3, ctx);
    ca_log(y, y, ctx);
    ca_div(x, x, y, ctx);
    ca_print(x, ctx);
    flint_printf("\n");

    flint_printf("6:  sqrt(2)*sqrt(3) - sqrt(6)                    = ");
    ca_sqrt_ui(x, 2, ctx);
    ca_sqrt_ui(y, 3, ctx);
    ca_sqrt_ui(z, 6, ctx);
    ca_mul(x, x, y, ctx);
    ca_sub(x, x, z, ctx);
    ca_print(x, ctx);
    flint_printf("\n");

    flint_printf("7:  exp(1+sqrt(2))*exp(1-sqrt(2))/(exp(1)^2)     = ");
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
    ca_print(x, ctx);
    flint_printf("\n");

    ca_clear(x, ctx);
    ca_clear(y, ctx);
    ca_clear(z, ctx);
    ca_ctx_clear(ctx);

    flint_printf("\n");
    TIMEIT_ONCE_STOP
    SHOW_MEMORY_USAGE

    flint_cleanup();
    return EXIT_SUCCESS;
}
