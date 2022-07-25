
#include "acb_theta.h"


int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("eld_points....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000 * arb_test_multiplier(); iter++)
      {
	slong g = 1 + n_randint(state, 4);
	slong d = 1 + n_randint(state, g);
	acb_theta_eld_t E;
	arb_mat_t Y;
	arb_t normsqr;
	arb_ptr offset;
	slong* last_coords;
	ulong a = n_randint(state, n_pow(2, g));
	ulong a_shift;
	slong prec = ACB_THETA_ELD_DEFAULT_PREC;
	slong mag_bits = n_randint(state, 2);
	slong k, j;
	slong try;
	slong* all_pts;
	slong* pt;
	int res;
	arb_mat_t vec;
	arb_t sqr, sum;
	
	acb_theta_eld_init(E, d, g);
	arb_mat_init(Y, g, g);
	arb_init(normsqr);
	offset = _arb_vec_init(g);
	last_coords = flint_malloc((g-d) * sizeof(slong));
	pt = flint_malloc(g * sizeof(slong));
	arb_mat_init(vec, g, 1);
	arb_init(sqr);
	arb_init(sum);

	arb_mat_randtest_cho(Y, state, prec, mag_bits);
	arb_randtest_pos(normsqr, state, prec, mag_bits);
	arb_mul_si(normsqr, normsqr, 1 + n_randint(state, 10), prec);

	a_shift = a;
	for (k = g-d-1; k >= 0; k--)
	  {	    
	    last_coords[k] = 2*n_randint(state, 5) + (a_shift % 2);
	    a_shift = a_shift >> 1;
	  }
	for (k = 0; k < g; k++) arb_randtest_precise(&offset[k], state, prec, mag_bits);
	
	acb_theta_eld_fill(E, Y, normsqr, offset, last_coords, a, prec);
	all_pts = flint_malloc(acb_theta_eld_nb_pts(E) * g * sizeof(slong));
	acb_theta_eld_points(all_pts, E);

	/* Test:
	   - all ellipsoid points must be within the box
	   - all ellipsoid points must have correct last coordinates
	   Then, generate random points:
	   - points inside ellipsoid must appear in all_pts
	   - points outside ellipsoid must have norm greater than normsqr
	*/

	for (k = 0; k < acb_theta_eld_nb_pts(E); k++)
	  {
	    for (j = 0; j < d; j++)
	      {
		if (FLINT_ABS(all_pts[k*g+j]) > acb_theta_eld_box(E, j))
		  {
		    flint_printf("FAIL: point outside box\n");
		    for (j = 0; j < g; j++) flint_printf("%wd ", all_pts[k*g+j]);
		    flint_printf("\nBox:\n");
		    for (j = 0; j < g; j++) flint_printf("%wd ", acb_theta_eld_box(E,j));
		    flint_printf("\n");
		    fflush(stdout);
		    flint_abort();
		  }
	      }
	    for (j = d; j < g; j++)
	      {
		if (all_pts[k*g+j] != acb_theta_eld_coord(E, j))
		  {
		    flint_printf("FAIL: incorrect coordinate\n");
		    for (j = 0; j < g; j++) flint_printf("%wd ", pt[j]);
		    fflush(stdout);
		    flint_abort();
		  }
	      }
	  }
	
	for (try = 0; try < 100; try++)
	  {
	    a_shift = a;
	    for (k = g-1; k >= 0; k--)
	      {		
		if (k >= d) pt[k] = last_coords[k-d];
		else
		  {
		    pt[k] = 2*n_randint(state, 2 + acb_theta_eld_box(E, k)/2);
		    pt[k] += (a_shift % 2);
		  }
		a_shift = a_shift >> 1;
	      }
	    if (acb_theta_eld_contains(E, pt))
	      {
		for (k = 0; k < acb_theta_eld_nb_pts(E); k++)
		  {
		    res = 1;
		    for (j = 0; j < g; j++)
		      {
			if (all_pts[k*g+j] != pt[j]) {res = 0; break;}
		      }
		    if (res == 1) break;
		  }
		if (!res)
		  {
		    flint_printf("FAIL: point not listed:\n");
		    for (j = 0; j < g; j++) flint_printf("%wd ", pt[j]);
		    fflush(stdout);
		    flint_abort();
		  }
	      }
	    
	    if (!acb_theta_eld_contains(E, pt))
	      {
		arb_mat_zero(vec);
		for (k = 0; k < d; k++) arb_set_si(arb_mat_entry(vec, k, 0), pt[k]);
		arb_mat_mul(vec, Y, vec, prec);
		arb_zero(sum);
		for (k = 0; k < d; k++)
		  {
		    arb_add(arb_mat_entry(vec, k, 0),
			    arb_mat_entry(vec, k, 0), &offset[k], prec);
		    arb_sqr(sqr, arb_mat_entry(vec, k, 0), prec);
		    arb_add(sum, sum, sqr, prec);
		  }
		if (arb_lt(sum, normsqr))
		  {
		    flint_printf("FAIL: small point not in ellipsoid\n");
		    for (j = 0; j < g; j++) flint_printf("%wd ", pt[j]);
		    flint_printf("\nCholesky:\n");
		    arb_mat_printd(Y, 10);
		    flint_printf("Norm of point: "); arb_printd(sum, 10);
		    flint_printf("\nCoordinates:\n");
		    for (j = 0; j < g; j++)
		      {
			arb_printd(arb_mat_entry(vec, j, 0), 10); flint_printf("\n");
		      }
		    flint_printf("Upper bound: "); arb_printd(normsqr, 10);
		    flint_printf("\na = %wu; total nb of points = %wd\n", a,
				 acb_theta_eld_nb_pts(E));
		    flint_printf("Offset:\n");
		    for (j = 0; j < g; j++)
		      {
			arb_printd(&offset[j], 10); flint_printf("\n");
		      }
		    flint_printf("Points:\n");
		    for (k = 0; k < acb_theta_eld_nb_pts(E); k++)
		      {
			for (j = 0; j < g; j++) flint_printf("%wd ", all_pts[k*g+j]);
			flint_printf("\n");
		      }
		    fflush(stdout);
		    flint_abort();		    
		  }
	      }
	  }
	 
	acb_theta_eld_clear(E);
	arb_mat_clear(Y);
	arb_clear(normsqr);
	_arb_vec_clear(offset, g);
	flint_free(last_coords);
	flint_free(all_pts);
	flint_free(pt);
	arb_mat_clear(vec);
	arb_clear(sqr);
	arb_clear(sum);
      }
    
    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
