
#include "acb_theta.h"

int main()
{
  slong iter;
  
  flint_printf("eld_points....");
  fflush(stdout);

  for (iter = 0; iter < 1; iter++)
    {
      arb_mat_t Y;
      arf_t epsilon;
      arb_t normsqr;
      arb_ptr offset;
      arb_eld_t E;
      slong* points;
      slong g = 2;
      ulong a = 0;
      slong prec = 40;
      slong nb;
      int v = 1;
      slong k, j;

      arb_mat_init(Y, g, g);
      arf_init(epsilon);
      arb_init(normsqr);
      offset = _arb_vec_init(g);
      arb_eld_init(E, g, g);

      arb_set_str(arb_mat_entry(Y,0,0), "0.9510565162", prec);
      arb_set_str(arb_mat_entry(Y,0,1), "0.363271264", prec);
      arb_set(arb_mat_entry(Y,1,0), arb_mat_entry(Y,0,1));
      arb_set_str(arb_mat_entry(Y,1,1), "0.9510565162", prec);
      arb_mat_cho(Y, Y, prec);
      arb_mat_transpose(Y, Y);
      arb_mat_scalar_mul_2exp_si(Y, Y, -1);
      
      if (v)
	{
	  flint_printf("Cholesky:\n");
	  arb_mat_printd(Y, 10);
	}
      
      arf_set_d(epsilon, 0.001);
      acb_theta_naive_radius(arb_midref(normsqr), Y, 0, epsilon, prec);
      _arb_vec_zero(offset, g);
      
      if (v)
	{
	  flint_printf("Predicted radius: ");
	  arb_printd(normsqr, 10); flint_printf("\n");
	  arf_set_d(arb_midref(normsqr), 7.3);
	  flint_printf("Enumeration radius: ");
	  arb_printd(normsqr, 10); flint_printf("\n");
	}

      arb_eld_fill(E, Y, normsqr, offset, NULL, a, prec);
      nb = arb_eld_nb_pts(E);
      flint_printf("Number of points: %wd\n", nb);
      points = flint_malloc(nb * g * sizeof(slong));
      arb_eld_points(points, E);
      for (k = 0; k < nb; k++)
	{
	  flint_printf("(%wd", points[k*g]);
	  for (j = 1; j < g; j++) flint_printf(",%wd", points[k*g+j]);
	  flint_printf(")\n");
	}

      arb_mat_clear(Y);
      arf_clear(epsilon);
      arb_clear(normsqr);
      _arb_vec_clear(offset, g);
      arb_eld_clear(E);
      flint_free(points);
    }
  
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
