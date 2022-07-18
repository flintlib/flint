
#include "acb_theta.h"


int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("eld_points....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100 * arb_test_multiplier(); iter++)
      {
	slong g = 1 + n_randint(state, 4);
	slong d = 1 + n_randint(state, g);
	arb_eld_t E;
	arb_mat_t Y;
	arb_t normsqr;
	arb_ptr offset;
	slong* last_coords;
	ulong a = n_randint(state, n_pow(2, g));
	slong prec = ARB_ELD_DEFAULT_PREC;
	slong mag_bits = n_randint(state, 2);
	slong k, j;
	slong try;
	slong* all_pts;
	slong* pt;
	int res;
	arb_mat_t vec;
	arb_t sqr, sum;
	
	arb_eld_init(E, d, g);
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
	arb_mul_si(normsqr, normsqr, 1 + n_randint(state, 5), prec);
	
	for (k = 0; k < g-d; k++) last_coords[k] = n_randint(state, 10);
	for (k = 0; k < g; k++) arb_randtest_precise(&offset[k], state, prec, mag_bits);
	
	arb_eld_fill(E, Y, normsqr, offset, last_coords, a, prec);
	all_pts = flint_malloc(arb_eld_nb_pts(E) * g * sizeof(slong));
	arb_eld_points(all_pts, E);

	/* Test:
	   - all ellipsoid points must be within the box
	   Then, generate random points:
	   - points inside ellipsoid must appear in all_pts
	   - points outside ellipsoid must have norm greater than normsqr
	*/

	for (k = 0; k < arb_eld_nb_pts(E); k++)
	  {
	    res = 1;
	    for (j = 0; j < g; j++)
	      {
		if (FLINT_ABS(all_pts[k*g+j]) > arb_eld_box(E,j+1))
		  {
		    flint_printf("FAIL: point outside box\n");
		    for (j = 0; j < g; j++) flint_printf("%wd ", pt[j]);
		    fflush(stdout);
		    flint_abort();
		  }
	      }
	  }
	
	for (try = 0; try < 100 * arb_test_multiplier(); try++)
	  {
	    for (k = g-1; k >= 0; k--)
	      {
		if (k >= d) pt[k] = last_coords[k-d];
		else
		  {
		    pt[k] = 2*n_randint(state, 1 + arb_eld_box(E, k+1)/2);
		    pt[k] += a % 2;
		  }
		a = a>>1;
	      }
	    if (arb_eld_contains(E, pt))
	      {
		for (k = 0; k < arb_eld_nb_pts(E); k++)
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
	    
	    if (!arb_eld_contains(E, pt))
	      {
		for (k = 0; k < g; k++) arb_set_si(arb_mat_entry(vec, k, 0), pt[k]);
		arb_mat_mul(vec, Y, vec, prec);
		arb_zero(sum);
		for (k = 0; k < g; k++)
		  {
		    arb_sub(arb_mat_entry(vec, k, 0),
			    arb_mat_entry(vec, k, 0), &offset[k], prec);
		    arb_sqr(sqr, arb_mat_entry(vec, k, 0), prec);
		    arb_add(sum, sum, sqr, prec);
		  }
		if (arb_lt(sum, normsqr))
		  {
		    flint_printf("FAIL: small point not in ellipsoid\n");
		    for (j = 0; j < g; j++) flint_printf("%wd ", pt[j]);
		    flint_printf("\nCholesky:");
		    arb_mat_printd(Y, 10);
		    flint_printf("Norm of point: "); arb_printd(sum, 10);		    
		    flint_printf("\nUpper bound: "); arb_printd(normsqr, 10);
		    flint_printf("\na = %wu; nb of points = %wd\n", a, arb_eld_nb_pts(E));
		    flint_printf("Offset:\n");
		    for (j = 0; j < g; j++)
		      {
			arb_printd(&offset[j], 10); flint_printf("\n");
		      }
		    flint_printf("Points:\n");
		    for (k = 0; k < arb_eld_nb_pts(E); k++)
		      {
			for (j = 0; j < g; j++) flint_printf("%wd ", all_pts[k*g+j]);
			flint_printf("\n");
		      }
		    fflush(stdout);
		    flint_abort();		    
		  }
	      }
	  }
	 
	arb_eld_clear(E);
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
