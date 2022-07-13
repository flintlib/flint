
#include "acb_theta.h"

void arb_eld_points(slong* pts, const arb_eld_t E)
{
  slong d = arb_eld_dim(E);
  slong g = arb_eld_ambient_dim(E);
  slong nr = arb_eld_nr(E);
  slong nl = arb_eld_nl(E);
  slong k, j, i;

  if (d == 1)
    {
      i = 0;
      for (k = arb_eld_min(E); k <= arb_eld_max(E); k += arb_eld_step(E))
	{
	  pts[i] = k;
	  for (j = 1; j < g; j++)
	    {
	      pts[i + j] = arb_eld_coord(E, j);
	    }
	  i += g;
	}
    }
  else /* d > 1 */
    {
      i = 0;
      for (k = 0; k < nr; k++)
	{
	  arb_eld_points(&pts[i], arb_eld_rchild(E, k));
	  i += g * arb_eld_nb_pts(arb_eld_rchild(E, k));
	}
      for (k = 0; k < nl; k++)
	{
	  arb_eld_points(&pts[i], arb_eld_lchild(E, k));
	  i += g * arb_eld_nb_pts(arb_eld_lchild(E, k));
	}
    }
}
