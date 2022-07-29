
#include "acb_theta.h"

void acb_theta_bound(arf_t rad, arf_t bound, acb_srcptr z, const acb_mat_t tau, slong prec)
{
  slong g = acb_mat_nrows(tau);
  arb_mat_t im;
  acb_mat_t pert;
  arb_t lambda;
  arf_t up;
  arb_mat_t z_pert;
  arb_mat_t prod;
  slong j, k;
  int res;

  arb_mat_init(im, g, g);
  acb_mat_init(pert, g, g);
  arb_init(lambda);
  arf_init(up);
  arb_mat_init(z_pert, g, 1);
  arb_mat_init(prod, 1, 1);
  
  acb_mat_get_imag(im, tau);
  
  /* Get lower bound on radius around tau */
  arb_mat_pos_radius(rad, im, prec);
  arf_mul_2exp_si(rad, rad, -1);

  /* Get upper bound for exponential sum */
  acb_mat_set(pert, tau);
  for (j = 0; j < g; j++)
    {
      for (k = 0; k < g; k++) acb_add_error_arf(acb_mat_entry(pert,j,k), rad);
    }
  acb_mat_get_imag(im, pert);
  arb_mat_pos_lambda(lambda, im, prec);
  arb_sqrt(lambda, lambda, prec);
  arb_inv(lambda, lambda, prec);
  arb_add_si(lambda, lambda, 1, prec);
  arb_pow_ui(lambda, lambda, g, prec);
  arb_get_ubound_arf(bound, lambda, prec);

  /* Multiply by upper bound for exponential term */
  for (k = 0; k < g; k++)
    {
      arb_set(arb_mat_entry(z_pert, k, 0), acb_imagref(&z[k]));
      arb_add_error_arf(arb_mat_entry(z_pert, k, 0), rad);
    }
  res = arb_mat_inv(im, im, prec);
  if (!res) arf_pos_inf(bound);

  arb_mat_mul(z_pert, im, z_pert, prec);
  arb_mat_mul(prod, z_pert, z_pert, prec);
  arb_const_pi(lambda, prec);
  arb_mul(lambda, lambda, arb_mat_entry(prod, 0, 0), prec);
  arb_exp(lambda, lambda, prec);
  arb_get_ubound_arf(up, lambda, prec);
  arf_mul(bound, bound, up, prec, ARF_RND_CEIL);
  
  arb_mat_clear(im);
  acb_mat_clear(pert);
  arb_clear(lambda);
  arf_clear(up);
  arb_mat_clear(z_pert);
  arb_mat_clear(prod);  
}
