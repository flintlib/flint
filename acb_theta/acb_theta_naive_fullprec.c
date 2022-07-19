
#include "acb_theta.h"

slong acb_theta_naive_fullprec(const arb_eld_t E, slong prec)
{
  return prec + ceil(ACB_THETA_NAIVE_FULLPREC_ADDLOG * n_flog(1 + arb_eld_nb_pts(E), 2));
}
