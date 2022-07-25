
#include "acb_theta.h"

void acb_theta_eld_points(slong* pts, const acb_theta_eld_t E)
{
  slong d = acb_theta_eld_dim(E);
  slong g = acb_theta_eld_ambient_dim(E);
  slong nr = acb_theta_eld_nr(E);
  slong nl = acb_theta_eld_nl(E);
  slong k, j, i;

  if (d == 1)
    {
      i = 0;
      for (k = acb_theta_eld_min(E); k <= acb_theta_eld_max(E); k += acb_theta_eld_step(E))
	{
	  pts[i] = k;
	  for (j = 1; j < g; j++)
	    {
	      pts[i + j] = acb_theta_eld_coord(E, j);
	    }
	  i += g;
	}
    }
  else /* d > 1 */
    {
      i = 0;
      for (k = 0; k < nr; k++)
	{
	  acb_theta_eld_points(&pts[i], acb_theta_eld_rchild(E, k));
	  i += g * acb_theta_eld_nb_pts(acb_theta_eld_rchild(E, k));
	}
      for (k = 0; k < nl; k++)
	{
	  acb_theta_eld_points(&pts[i], acb_theta_eld_lchild(E, k));
	  i += g * acb_theta_eld_nb_pts(acb_theta_eld_lchild(E, k));
	}
    }
}
