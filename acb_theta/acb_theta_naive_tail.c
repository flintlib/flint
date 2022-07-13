
#include "acb_theta.h"

/* Cf note */

void acb_theta_naive_tail(arf_t B, const arf_t R, const arb_mat_t Y, slong p, slong prec)
{
  arb_t res, temp;
  arb_t Rmod;
  slong g = arb_mat_nrows(Y);
  slong k;

  arb_init(res);
  arb_init(temp);
  arb_init(Rmod);
  
  /* Ensure assumptions R\geq 4, R\geq 2p are satisfied */
  arb_set_arf(Rmod, R);
  arb_set_si(temp, FLINT_MAX(4, 2*p));
  arb_max(Rmod, Rmod, temp, prec);
  
  /* Evaluate upper bound on tail */
  arb_one(res);
  arb_mul_2exp_si(res, res, 2*g+2);

  arb_sqrt(temp, Rmod, prec);
  arb_pow_ui(temp, temp, g-1+2*p, prec);
  arb_mul(res, res, temp, prec);

  arb_neg(temp, Rmod);
  arb_exp(temp, temp, prec);
  arb_mul(res, res, temp, prec);

  for (k = 0; k < g; k++)
    {
      arb_inv(temp, arb_mat_entry(Y, k, k), prec);
      arb_add_si(temp, temp, 1, prec);
      arb_mul(res, res, temp, prec);
    }
  arb_get_ubound_arf(B, res, prec);
  
  arb_clear(res);
  arb_clear(temp);
  arb_clear(Rmod);
}
