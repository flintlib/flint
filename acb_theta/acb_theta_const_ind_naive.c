
#include "acb_theta.h"

static void set_precomp(arb_eld_t E, acb_theta_precomp_t D, arf_t epsilon,
			ulong ab, const acb_mat_t tau, slong prec)
{
  arf_t R;
  arb_mat_t Y;
  arb_t pi;
  arb_t normsqr;
  arb_ptr offset;
  slong g = acb_mat_nrows(tau);
  slong eld_prec = ARB_ELD_DEFAULT_PREC;
  int res;
  slong k;

  arf_init(R);
  arb_mat_init(Y, g, g);
  arb_init(normsqr);
  arb_init(pi);
  offset = _arb_vec_init(g);
  
  arf_one(epsilon);
  arf_mul_2exp_si(epsilon, epsilon, -prec + ACB_THETA_NAIVE_EPS_2EXP);

  acb_mat_get_imag(Y, tau);
  arb_const_pi(pi, prec);
  arb_mat_scalar_mul_arb(Y, Y, pi, prec);
  
  res = arb_mat_cho(Y, Y, eld_prec);
  if (!res)
    {
      eld_prec = prec;
      arb_mat_cho(Y, Y, eld_prec);
    }
  if (!res)
    {
      flint_printf("Error: not positive definite\n");
      fflush(stdout);
      flint_abort();
    }
  arb_mat_transpose(Y, Y);
  acb_theta_naive_radius(R, Y, 0, epsilon, eld_prec);

  flint_printf("Cholesky diagonal:\n");
  for (k = 0; k < g; k++)
    {
      arb_printd(arb_mat_entry(Y,k,k), 10); flint_printf("\n");
    }

  flint_printf("Chosen error and ellipsoid radius:\n");
  arf_printd(epsilon, 10); flint_printf("\n");
  arf_printd(R, 10); flint_printf("\n");
  
  arb_set_arf(normsqr, R);
  arb_mul_2exp_si(normsqr, normsqr, 2);
  
  _arb_vec_zero(offset, g);
  
  arb_eld_fill(E, Y, normsqr, offset, NULL, ab >> g, eld_prec);  
  acb_theta_precomp_set(D, tau, E, prec);

  flint_printf("Ellipsoid box:");
  for (k = 1; k <= g; k++) flint_printf(" %wd", arb_eld_box(E, k));
  flint_printf("\n");

  arf_clear(R);
  arb_mat_clear(Y);
  arb_clear(normsqr);
  _arb_vec_clear(offset, g);  
}

static void worker_dim0(acb_ptr th, const acb_t term, slong* coords, slong g,
			ulong ab, slong ord, slong prec, slong fullprec)
{
  acb_t x;
  slong sgn = 0;
  slong k;

  acb_init(x);

  
  flint_printf("Exponential term");
  for (k = 0; k < g; k++) flint_printf(" %wd", coords[k]);
  flint_printf(":\n");
  acb_printd(term, 30);
  flint_printf("\n");
  
  for (k = 0; k < g; k++)
    {
      if (ab & 1)
	{
	  sgn += 4 + coords[g-1-k] % 4; /* & 3 ? */
	}
      ab = ab >> 1;
    }
  sgn = sgn % 4;
  
  acb_set(x, term);
  if (sgn == 1) acb_mul_onei(x, x);
  else if (sgn == 2) acb_neg(x, x);
  else if (sgn == 3) acb_div_onei(x, x);

  acb_add(th, th, x, fullprec);    
  acb_clear(x);
}

void acb_theta_const_ind_naive(acb_t th, ulong ab, const acb_mat_t tau, slong prec)
{
  arb_eld_t E;
  acb_theta_precomp_t D;
  arf_t epsilon;
  acb_mat_t lin_powers;
  acb_ptr exp_z;
  acb_t cofactor;
  slong ord = 0;
  slong fullprec;
  slong g = acb_mat_nrows(tau);
  slong k;

  arb_eld_init(E, g, g);
  acb_theta_precomp_init(D, g);
  arf_init(epsilon);
  acb_mat_init(lin_powers, g, g);
  exp_z = _acb_vec_init(g);
  acb_init(cofactor);

  set_precomp(E, D, epsilon, ab, tau, prec);

  acb_mat_set(lin_powers, acb_theta_precomp_exp_mat(D));
  for (k = 0; k < g; k++) acb_one(&exp_z[k]);
  acb_one(cofactor);
  fullprec = prec + ceil(ACB_THETA_NAIVE_FULLPREC_ADDLOG * n_flog(1 + arb_eld_nb_pts(E), 2));

  acb_zero(th);
  acb_theta_naive_worker_rec(th, lin_powers, E, D, exp_z, cofactor, ab, ord,
			     prec, fullprec, worker_dim0);

  arb_eld_clear(E);
  acb_theta_precomp_clear(D);
  arf_clear(epsilon);
  acb_mat_clear(lin_powers);
  _acb_vec_clear(exp_z, g);
  acb_clear(cofactor);  
}
