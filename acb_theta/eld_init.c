
#include "acb_theta.h"

void acb_theta_eld_init(acb_theta_eld_t E, slong d, slong g)
{
  acb_theta_eld_dim(E) = d;
  acb_theta_eld_ambient_dim(E) = g;
  E->last_coords = flint_malloc((g-d) * sizeof(slong));
  acb_theta_eld_offset(E) = _arb_vec_init(d);
  arb_init(acb_theta_eld_normsqr(E));

  arb_init(acb_theta_eld_ctr(E));
  arb_init(acb_theta_eld_rad(E));
  E->rchildren = NULL;
  acb_theta_eld_nr(E) = 0;
  E->lchildren = NULL;
  acb_theta_eld_nl(E) = 0;
  E->box = flint_malloc(d * sizeof(slong));
}

