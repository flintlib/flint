
#include "acb_theta.h"

void acb_theta_agm_collect(acb_ptr all_r0, arb_t M0, arb_t minf, arb_ptr mi,
			   slong nb_bad, const acb_mat_t tau, const fmpz_mat_t N, slong prec)
{
  acb_mat_t z;
  arb_t abs;
  slong g = acb_mat_nrows(tau);
  slong n = 1<<g;
  slong k, j;
  slong lowprec = ACB_THETA_AGM_LOWPREC;
  
  acb_mat_init(z, g, g);
  arb_init(abs);
  
  acb_siegel_transform(z, N, tau, prec);

  for (k = 0; k < nb_bad; k++)
    {
      acb_theta_naive_const(all_r0 + k*n, z, lowprec);
      acb_mat_scalar_mul_2exp_si(z, z, 1);
    }
  
  arb_zero(M0);
  for (k = 0; k < n*nb_bad; k++)
    {
      arb_max(M0, M0, &all_r0[k], lowprec);
    }
  arb_sqr(M0, M0, lowprec);

  arb_set_si(minf, 19);
  arb_div_si(minf, minf, 20, lowprec); /* Cf. agm_nb_bad_steps */

  for (k = 0; k < nb_bad; k++)
    {
      arb_pos_inf(&mi[k]);
      for (j = 0; j < n; j++)
	{
	  acb_abs(abs, &all_r0[k*n + j], lowprec);	  
	  arb_min(&mi[k], &mi[k], abs, lowprec);
	}
    }
  
  acb_mat_clear(z);
  arb_clear(abs);  
}
