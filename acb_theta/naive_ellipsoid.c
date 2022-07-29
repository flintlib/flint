
#include "acb_theta.h"

void acb_theta_naive_ellipsoid(acb_theta_eld_t E, arf_t epsilon,
			       ulong ab, int all, int unif, slong ord,
			       acb_srcptr z, const acb_mat_t tau, slong prec)
{  
  arf_t R;
  arb_mat_t im;
  arb_mat_t cho;
  arb_t pi;
  arb_t normsqr;
  arb_mat_t imz;
  slong* translate;
  arb_ptr offset;
  slong g = acb_mat_nrows(tau);
  slong eld_prec = ACB_THETA_ELD_DEFAULT_PREC;
  int res;
  slong k;

  arf_init(R);
  arb_mat_init(im, g, g);
  arb_mat_init(cho, g, g);
  arb_init(normsqr);
  arb_init(pi);
  arb_mat_init(imz, g, 1);
  offset = _arb_vec_init(g);
  translate = flint_malloc(g * sizeof(slong));
  
  arf_one(epsilon);
  arf_mul_2exp_si(epsilon, epsilon, -prec + ACB_THETA_NAIVE_EPS_2EXP);

  acb_mat_get_imag(im, tau);
  arb_const_pi(pi, prec);
  arb_mat_scalar_mul_arb(cho, im, pi, prec);
  
  res = arb_mat_cho(cho, cho, eld_prec);
  if (!res)
    {
      eld_prec = prec;
      arb_mat_cho(cho, cho, eld_prec);
    }
  if (!res)
    {
      flint_printf("acb_theta_naive_ellipsoid: imaginary part is not positive definite\n");
      fflush(stdout);
      flint_abort();
    }  
  arb_mat_transpose(cho, cho);
  
  if (all) /* need all points in Z^g */
    {
      ab = 0;
      arb_mat_scalar_mul_2exp_si(cho, cho, -1);
    }  
  acb_theta_naive_radius(R, cho, ord, epsilon, eld_prec);
  
  arb_set_arf(normsqr, R);
  arb_mul_2exp_si(normsqr, normsqr, 2);
  
  if (unif) /* any offset less than 1/2 */
    {
      arb_one(pi);
      arb_mul_2exp_si(pi, pi, -1);
      for (k = 0; k < g; k++)
	{
	  arb_unit_interval(&offset[k]);
	  arb_sub(&offset[k], &offset[k], pi, eld_prec);
	}
    }
  else /* set offset in terms of z */
    {
      for (k = 0; k < g; k++)
	{
	  arb_set(arb_mat_entry(imz, k, 0), acb_imagref(&z[k]));	  
	}
      arb_mat_inv(im, im, eld_prec);
      arb_mat_mul(imz, im, imz, prec);
      arb_mat_mul(imz, cho, imz, prec);
      for (k = 0; k < g; k++)
	{
	  arb_set(&offset[k], arb_mat_entry(imz, k, 0));
	}
    }
  
  acb_theta_eld_fill(E, cho, normsqr, offset, NULL, ab >> g, eld_prec);
  
  /* exponential error factor in terms of z */
  for (k = 0; k < g; k++)
    {
      arb_set(arb_mat_entry(imz, k, 0), acb_imagref(&z[k]));	  
    }
  arb_mat_mul(imz, im, imz, prec);
  arb_zero(normsqr);
  for (k = 0; k < g; k++)
    {
      arb_sqr(pi, arb_mat_entry(imz, k, 0), prec);
      arb_add(normsqr, normsqr, pi, prec);
    }
  arb_const_pi(pi, prec);
  arb_mul(normsqr, normsqr, pi, prec);
  arb_exp(normsqr, normsqr, prec);
  arb_get_ubound_arf(R, normsqr, prec);
  arf_mul(epsilon, epsilon, R, prec, ARF_RND_CEIL);
    
  arf_clear(R);
  arb_mat_clear(im);
  arb_mat_clear(cho);  
  arb_clear(normsqr);
  arb_clear(pi);
  arb_mat_clear(imz);
  _arb_vec_clear(offset, g);
  flint_free(translate);
}
