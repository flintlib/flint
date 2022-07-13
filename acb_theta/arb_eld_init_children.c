
#include "acb_theta.h"

void arb_eld_init_children(arb_eld_t E, slong nr, slong nl)
{
  slong d = arb_eld_dim(E);
  slong g = arb_eld_ambient_dim(E);
  slong k;

  if (nr > 0)
    {
      E->rchildren = flint_malloc(nr * sizeof(struct arb_eld_struct));
      arb_eld_nr(E) = nr;
      for (k = 0; k < nr; k++) arb_eld_init(arb_eld_rchild(E, k), d-1, g);
    }
  if (nl > 0)
    {
      E->lchildren = flint_malloc(nl * sizeof(struct arb_eld_struct));
      arb_eld_nl(E) = nl;
      for (k = 0; k < nl; k++) arb_eld_init(arb_eld_lchild(E, k), d-1, g);
    }
}
