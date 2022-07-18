
#include "acb_theta.h"

slong acb_theta_naive_newprec(slong prec, slong coord, slong dist, slong max_dist,
			      slong step, slong ord)
{
  double r = ((double) dist)/(max_dist + step);
  double neg = ACB_THETA_NAIVE_NEWPREC_MARGIN * r * r * prec;
  double pos = ord * n_clog(1 + FLINT_ABS(coord), 2);

  return ceil((double) prec - neg + pos);
}
