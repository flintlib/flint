
#include "acb_theta.h"

void acb_theta_eld_init_children(acb_theta_eld_t E, slong nr, slong nl)
{
  slong d = acb_theta_eld_dim(E);
  slong g = acb_theta_eld_ambient_dim(E);
  slong k;

  if (nr > 0)
    {
      E->rchildren = flint_malloc(nr * sizeof(struct acb_theta_eld_struct));
      acb_theta_eld_nr(E) = nr;
      for (k = 0; k < nr; k++) acb_theta_eld_init(acb_theta_eld_rchild(E, k), d-1, g);
    }
  if (nl > 0)
    {
      E->lchildren = flint_malloc(nl * sizeof(struct acb_theta_eld_struct));
      acb_theta_eld_nl(E) = nl;
      for (k = 0; k < nl; k++) acb_theta_eld_init(acb_theta_eld_lchild(E, k), d-1, g);
    }
}
