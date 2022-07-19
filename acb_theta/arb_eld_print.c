
#include "acb_theta.h"

void arb_eld_print(const arb_eld_t E)
{
  slong d = arb_eld_dim(E);
  slong g = arb_eld_ambient_dim(E);
  slong k;

  for (k = 0; k < g-d; k++) flint_printf("  ");
  flint_printf("Slice (...");
  for (k = 0; k < g-d; k++) flint_printf(", %wd", arb_eld_coord(E, k+d));
  flint_printf("): from %wd to %wd by %wd (mid: %wd)\n",
	       arb_eld_min(E),
	       arb_eld_max(E),
	       arb_eld_step(E),
	       arb_eld_mid(E));
  if (d > 1)
    {
      for (k = 0; k < arb_eld_nr(E); k++)
	{
	  arb_eld_print(arb_eld_rchild(E,k));
	}
      for (k = 0; k < arb_eld_nl(E); k++)
	{
	  arb_eld_print(arb_eld_lchild(E,k));
	}
    }
}
