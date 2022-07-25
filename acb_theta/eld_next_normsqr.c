
#include "acb_theta.h"

void acb_theta_eld_next_normsqr(arb_t next_normsqr, const arb_t normsqr, const arb_t gamma,
				const arb_t ctr, slong k, slong prec)
{
  arb_t x;
  arb_init(x);

  /* Set next_normsqr to normsqr - gamma^2(ctr - k)^2 */
  arb_sub_si(x, ctr, k, prec);
  arb_mul(x, x, gamma, prec);
  arb_sqr(x, x, prec);
  arb_sub(next_normsqr, normsqr, x, prec);

  arb_clear(x);
}
