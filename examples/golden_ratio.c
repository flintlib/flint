/* This file is public domain. Author: Raoul Bourquin. */

#include <stdlib.h>
#include "ca.h"


void main_fexpr()
{
    fexpr_t Phi;
    fexpr_init(Phi);

    flint_printf("Evaluating Phi as fexpr:\n");

    fexpr_set_symbol_str(Phi, "GoldenRatio");

    fexpr_print(Phi);
    printf("\n\n");

    fexpr_clear(Phi);
}


void main_ca()
{
    ca_ctx_t ctx;
    ca_t Phi;
    ca_ctx_init(ctx);
    ca_init(Phi, ctx);

    flint_printf("Evaluating Phi as ca:\n");

    ca_phi(Phi, ctx);

    ca_print(Phi, ctx);
    printf("\n\n");

    ca_clear(Phi, ctx);
}


void main_qqbar()
{
    qqbar_t Phi;
    qqbar_init(Phi);

    flint_printf("Evaluating Phi as qqbar:\n");

    qqbar_phi(Phi);

    qqbar_printn(Phi, 50);
    printf("\n");

    qqbar_clear(Phi);
}


int main(int argc, char *argv[])
{
  main_fexpr();
  main_ca();
  main_qqbar();

  flint_cleanup();
  return EXIT_SUCCESS;
}
