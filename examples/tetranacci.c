/* This file is public domain. Author: Raoul Bourquin. */

#include <stdlib.h>
#include "ca.h"


void main_fexpr()
{
    fexpr_t T;
    fexpr_init(T);

    flint_printf("Evaluating Tt as fexpr:\n");

    fexpr_set_symbol_str(T, "TetranacciConstant");

    fexpr_print(T);
    printf("\n\n");

    fexpr_clear(T);
}


void main_ca()
{
    ca_ctx_t ctx;
    ca_t T;
    ca_ctx_init(ctx);
    ca_init(T, ctx);

    flint_printf("Evaluating Tt as ca:\n");

    ca_tetranacci_constant(T, ctx);

    ca_print(T, ctx);
    printf("\n\n");

    ca_clear(T, ctx);
}


void main_qqbar()
{
    qqbar_t T;
    qqbar_init(T);

    flint_printf("Evaluating Tt as qqbar:\n");

    qqbar_tetranacci_constant(T);

    qqbar_printn(T, 50);
    printf("\n");

    qqbar_clear(T);
}


int main(int argc, char *argv[])
{
  main_fexpr();
  main_ca();
  main_qqbar();

  flint_cleanup();
  return EXIT_SUCCESS;
}
