
#include "acb_theta.h"

void acb_theta_balance(acb_ptr z2, acb_mat_t tau2, fmpz_mat_t mat,
	acb_srcptr z, const acb_mat_t tau, slong j)
{
    slong g = acb_mat_nrows(tau);    
    slong k, l;

    fmpz_mat_one(mat);
    
    for (k = 0; k <= j; k++)
    {
	fmpz_zero(fmpz_mat_entry(mat, k, k));
	fmpz_zero(fmpz_mat_entry(mat, k+g, k+g));
	fmpz_one(fmpz_mat_entry(mat, k+g, k));
	fmpz_set_si(fmpz_mat_entry(mat, k, k+g), -1);
    }
    
    acb_mat_set(tau2, tau);
    
    for (k = 0; k <= j; k++)
    {
	for (l = 0; l <= j; l++)
	{
	    acb_mul_2exp_si(acb_mat_entry(tau2, k, l),
		    acb_mat_entry(tau2, k, l), 1);
	}
    }
    for (k = j+1; k < g; k++)
    {
	for (l = j+1; l < g; l++)
	{
	    acb_mul_2exp_si(acb_mat_entry(tau2, k, l),
		    acb_mat_entry(tau2, k, l), -1);
	}
    }

    _acb_vec_set(z2, z, g);
    for (k = 0; k <= j; k++)
    {
	acb_mul_2exp_si(&z2[k], &z2[k], 1);
    }
}
