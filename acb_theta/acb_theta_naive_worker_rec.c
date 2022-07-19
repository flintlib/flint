
#include "acb_theta.h"

void acb_theta_naive_worker_rec(acb_ptr th, acb_mat_t lin_powers,
				const arb_eld_t E, const acb_theta_precomp_t D,
				acb_srcptr exp_z, const acb_t cofactor,
				ulong ab, slong ord, slong prec, slong fullprec,
				acb_theta_naive_worker_t worker_dim0)
{  
  slong d = arb_eld_dim(E);
  slong g = arb_eld_ambient_dim(E);
  slong nr = arb_eld_nr(E);
  slong nl = arb_eld_nl(E);
  slong min = arb_eld_min(E);
  slong mid = arb_eld_mid(E);
  slong max = arb_eld_max(E);
  slong step = arb_eld_step(E);
  acb_t start_cf, diff_cf, lin_cf, full_cf; /* Set up next cofactor */
  acb_ptr start_lin_powers, diff_lin_powers; /* Set up next lin_powers */
  slong newprec;
  slong k, j, c;

  /* Catch cases: no points in ellipsoid; d=1 */
  if (arb_eld_nb_pts(E) == 0)
    {
      return;
    }
  else if (d == 1)
    {
      acb_init(lin_cf);
      acb_set(lin_cf, &exp_z[0]);
      for (k = 1; k < g; k++)
	{
	  acb_mul(lin_cf, lin_cf,
		  acb_mat_entry(lin_powers, 0, k), prec);
	}      
      acb_theta_naive_worker_dim1(th, E, D, lin_cf, cofactor,
				  ab, ord, prec, fullprec, worker_dim0);
      acb_clear(lin_cf);
      return;
    }

  acb_init(start_cf);
  acb_init(diff_cf);
  acb_init(lin_cf);
  acb_init(full_cf);
  start_lin_powers = _acb_vec_init(d-1);
  diff_lin_powers = _acb_vec_init(d-1);
  
  /* Set up things for new cofactor */
  acb_set(diff_cf, &exp_z[d-1]);
  for (k = d; k < g; k++)
    {
      acb_mul(diff_cf, diff_cf, acb_mat_entry(lin_powers, d-1, k), prec);
    }
  acb_pow_si(start_cf, diff_cf, mid, prec);
  acb_mul(start_cf, start_cf, cofactor, prec);
  acb_pow_si(diff_cf, diff_cf, step, prec);

  /* Set up things to update entries (k,d) of lin_powers, k < d */
  for (k = 0; k < d-1; k++)
    {
      acb_pow_si(&diff_lin_powers[k],
		 acb_mat_entry(acb_theta_precomp_exp_mat(D), k, d-1), step, prec);
      acb_pow_si(&start_lin_powers[k],
		 acb_mat_entry(acb_theta_precomp_exp_mat(D), k, d-1), mid, prec);
    }

  /* Right loop */
  acb_set(lin_cf, start_cf);
  for (k = 0; k < d-1; k++)
    {
      acb_set(acb_mat_entry(lin_powers, k, d-1), &start_lin_powers[k]);
    }  
  for (k = 0; k < nr; k++)
    {
      c = mid + k*step;
      newprec = acb_theta_naive_newprec(prec, c, c-mid, max-mid, step, ord);
      if (k > 0) /* Update lin_cf, lin_powers using diff */
	{
	  for (j = 0; j < d-1; j++)
	    {
	      acb_mul(acb_mat_entry(lin_powers, j, d-1),
		      acb_mat_entry(lin_powers, j, d-1), &diff_lin_powers[j], newprec);
	    }
	  acb_mul(lin_cf, lin_cf, diff_cf, newprec);
	}
      
      acb_mul(full_cf, lin_cf, acb_theta_precomp_sqr_pow(D, d, FLINT_ABS(c)/step), newprec);
      acb_theta_naive_worker_rec(th, lin_powers, arb_eld_rchild(E,k), D, exp_z, full_cf,
				 ab, ord, newprec, fullprec, worker_dim0);
    }

  /* Left loop */
  acb_set(lin_cf, start_cf);
  for (k = 0; k < d-1; k++)
    {
      acb_set(acb_mat_entry(lin_powers, k, d-1), &start_lin_powers[k]);
    }
  acb_inv(diff_cf, diff_cf, prec);
  for (k = 0; k < d-1; k++)
    {
      acb_inv(&diff_lin_powers[k], &diff_lin_powers[k], prec);
    }
  for (k = 0; k < nl; k++)
    {      
      c = mid - (k+1)*step;
      newprec = acb_theta_naive_newprec(prec, c, mid-c, mid-min, step, ord);
      for (j = 0; j < d-1; j++)
	{
	  acb_mul(acb_mat_entry(lin_powers, j, d-1),
		  acb_mat_entry(lin_powers, j, d-1), &diff_lin_powers[j], newprec);
	}
      acb_mul(lin_cf, lin_cf, diff_cf, newprec);
      
      acb_mul(full_cf, lin_cf, acb_theta_precomp_sqr_pow(D, d, FLINT_ABS(c)/step), newprec);
      acb_theta_naive_worker_rec(th, lin_powers, arb_eld_rchild(E,k), D, exp_z, full_cf,
				 ab, ord, newprec, fullprec, worker_dim0);
    }
  
  acb_clear(start_cf);
  acb_clear(diff_cf);
  acb_clear(lin_cf);
  acb_clear(full_cf);
  _acb_vec_clear(start_lin_powers, d-1);
  _acb_vec_clear(diff_lin_powers, d-1);
}
