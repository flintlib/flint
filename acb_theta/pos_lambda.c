
#include "acb_theta.h"

void
arb_mat_pos_lambda(arb_t lambda, const arb_mat_t mat, slong prec)
{
    arb_poly_t poly;

    arb_poly_init(poly);
  
    arb_mat_charpoly(poly, mat, prec);
    arb_div(lambda, arb_poly_get_coeff_ptr(poly, 0),
	    arb_poly_get_coeff_ptr(poly, 1), prec);
    arb_neg(lambda, lambda);

    arb_poly_clear(poly);
}
