
#include "acb_theta.h"

void arb_eld_clear(arb_eld_t E)
{
  slong k;
  slong nr = arb_eld_nr(E);
  slong nl = arb_eld_nl(E);
  
  if (nr > 0)
    {
      for (k = 0; k < nr; k++) arb_eld_clear(arb_eld_rchild(E, k));
      flint_free(E->rchildren);
    }
  if (nl > 0)
    {
      for (k = 0; k < nl; k++) arb_eld_clear(arb_eld_lchild(E, k));
      flint_free(E->lchildren);
    }

  flint_free(E->last_coords);
  _arb_vec_clear(arb_eld_offset(E), arb_eld_dim(E));
  arb_clear(arb_eld_normsqr(E));
  arb_clear(arb_eld_rad(E));
  arb_clear(arb_eld_ctr(E));
  flint_free(E->box);
}
