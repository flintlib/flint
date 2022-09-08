
#include "acb_theta.h"

void
acb_theta_precomp_set(acb_theta_precomp_t D, const acb_mat_t tau,
        const acb_theta_eld_t E, slong prec)
{
    slong g = acb_theta_eld_ambient_dim(E);
    arb_t pi4;
    acb_t c, dc, ddc;
    slong k, j, s;
    slong nb_pow;
    
    if (acb_theta_eld_nb_pts(E) == 0) return;
    
    arb_init(pi4);
    acb_init(c);
    acb_init(dc);
    acb_init(ddc);
    
    arb_const_pi(pi4, prec);
    arb_mul_2exp_si(pi4, pi4, -2);
    
    /* Set matrix of exponentials */
    for (k = 0; k < g; k++)
    {
        for (j = k; j < g; j++)
	{
            acb_mul_arb(c, acb_mat_entry(tau, k, j), pi4, prec);
            acb_mul_onei(c, c);
            if (k != j) acb_mul_2exp_si(c, c, 1);
            acb_exp(c, c, prec);
            acb_set(acb_mat_entry(acb_theta_precomp_exp_mat(D), k, j), c);
	}
    }
    
    /* Set indices */
    D->indices[0] = 0;
    for (k = 0; k < g; k++)
    {
        nb_pow = acb_theta_precomp_box(E, k) / 2 + 1;
        D->indices[k+1] = D->indices[k] + nb_pow;
    }

  /* Init and set square powers; addition chains unnecessary */
  D->sqr_powers = _acb_vec_init(D->indices[g]);
  for (k = 0; k < g; k++)
    {
      acb_set(ddc, acb_mat_entry(acb_theta_precomp_exp_mat(D), k, k));
      s = acb_theta_precomp_box(E, k) % 2;
      acb_pow_si(c, ddc, s, prec);
      acb_pow_si(dc, ddc, 4*s + 4, prec);
      acb_pow_si(ddc, ddc, 8, prec);
      for (j = 0; s + 2*j <= acb_theta_precomp_box(E, k); j++)
	{
	  acb_set(acb_theta_precomp_sqr_pow(D, k, j), c);
	  acb_mul(c, c, dc, prec);
	  acb_mul(dc, dc, ddc, prec);
	}
    }
  
  arb_clear(pi4);
  acb_clear(c);
  acb_clear(dc);
  acb_clear(ddc);
}
