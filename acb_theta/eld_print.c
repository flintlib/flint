
#include "acb_theta.h"

void
acb_theta_eld_print(const acb_theta_eld_t E)
{
  slong d = acb_theta_eld_dim(E);
  slong g = acb_theta_eld_ambient_dim(E);
  slong k;

  for (k = 0; k < g-d; k++) flint_printf("  ");
  flint_printf("Slice (...");
  for (k = 0; k < g-d; k++) flint_printf(", %wd", acb_theta_eld_coord(E, k+d));
  flint_printf("): from %wd to %wd by %wd (mid: %wd)\n",
	       acb_theta_eld_min(E),
	       acb_theta_eld_max(E),
	       2,
	       acb_theta_eld_mid(E));
  if (d > 1)
    {
      for (k = 0; k < acb_theta_eld_nr(E); k++)
	{
	  acb_theta_eld_print(acb_theta_eld_rchild(E,k));
	}
      for (k = 0; k < acb_theta_eld_nl(E); k++)
	{
	  acb_theta_eld_print(acb_theta_eld_lchild(E,k));
	}
    }
}
