
#include "acb_theta.h"

void arb_eld_init(arb_eld_t E, slong d, slong g)
{
  arb_eld_dim(E) = d;
  arb_eld_ambient_dim(E) = g;
  E->last_coords = flint_malloc((g-d) * sizeof(slong));
  arb_eld_offset(E) = _arb_vec_init(d);
  arb_init(arb_eld_normsqr(E));

  arb_init(arb_eld_ctr(E));
  arb_init(arb_eld_rad(E));
  E->rchildren = NULL;
  arb_eld_nr(E) = 0;
  E->lchildren = NULL;
  arb_eld_nl(E) = 0;
  E->box = flint_malloc(d * sizeof(slong));
}

