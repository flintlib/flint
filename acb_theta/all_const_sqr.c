
#include "acb_theta.h"

static void
fmpz_mat_balance(fmpz_mat_t D, slong j)
{
    slong g = acb_mat_nrows(D)/2;
    slong k;

    fmpz_mat_one(D);
    
    for (k = 0; k <= j; k++)
    {
	fmpz_zero(fmpz_mat_entry(D, k, k));
	fmpz_zero(fmpz_mat_entry(D, k+g, k+g));
	fmpz_one(fmpz_mat_entry(D, k+g, k));
	fmpz_set_si(fmpz_mat_entry(D, k, k+g), -1);
    }
}

static void
acb_mat_balance(acb_mat_t res, const acb_mat_t tau, slong j)
{
    slong g = acb_mat_nrows(tau);
    slong k, l;

    acb_mat_set(res, tau);
    
    for (k = 0; k <= j; k++)
    {
	for (l = 0; l <= j; l++)
	{
	    acb_mul_2exp_si(acb_mat_entry(res, k, l),
		    acb_mat_entry(res, k, l), 1);
	}
    }
    for (k = j+1; k < g; k++)
    {
	for (l = j+1; l < g; l++)
	{
	    acb_mul_2exp_si(acb_mat_entry(res, k, l),
		    acb_mat_entry(res, k, l), -1);
	}
    }
}

static int
is_balanced(acb_mat_t res, fmpz_mat_t D, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    arb_t test;
    slong j;
    int r = 1;

    arb_init(test);
    
    for (j = 0; j < g-1; j++)
    {
	arb_mul_2exp_si(test, acb_imagref(arb_mat_entry(tau, j, j)), 2);
	if (arb_lt(test, acb_imagref(arb_mat_entry(tau, j+1, j+1))))
	{
	    flint_printf("(is_balanced) Found unbalanced at j=%wd\n", j);
	    
	    r = 0;
	    fmpz_mat_balance(D, j);
	    acb_mat_balance(res, tau, j);
	    break;
	}
    }

    arb_clear(test);
    return r;
}

static int
accept_naive(const acb_mat_t tau, slong prec)
{
    arb_t test;
    int res;
    
    arb_init(test);

    arb_set(test, acb_imagref(acb_mat_entry(tau, 0, 0)));
    arb_mul_si(test, test, ACB_THETA_NAIVE_CMP, prec);
    arb_sub_si(test, test, prec, prec);

    res = arb_is_positive(test);

    arb_clear(test);
    return res;
}

static void
theta_naive(acb_ptr th2, const acb_mat_t tau, slong prec)
{    
    slong g = acb_mat_nrows(tau);
    slong k;
    
    acb_theta_naive_all_const(th2, tau, prec);
    for (k = 0; k < (1<<(2*g)); k++)
    {
	acb_sqr(&th2[k], &th2[k], prec);
    }
}

static int
lowprec_roots(acb_ptr roots, const acb_mat_t tau, const fmpz_mat_t mat,
	slong lowprec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1<<g;
    acb_ptr th;
    ulong k;
    int res = 1;

    th = _acb_vec_init(n*n);

    acb_theta_naive_all_const(th, tau, lowprec);
    acb_theta_transform_proj(roots, th, mat, lowprec);

    flint_printf("Lowprec roots:\n");
    for (k = 0; k < n; k++)
    {
	acb_printd(&roots[k], 10); flint_printf("\n");
    }

    for (k = 0; k < n; k++)
    {
	if (acb_contains_zero(&roots[k]))
	{
	    res = 0;
	    break;
	}
    }

    _acb_vec_clear(th, n*n);
    return res;
}

void
acb_theta_all_const_sqr(acb_ptr th2, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1<<g;
    fmpz_mat_t D, R;
    acb_mat_t aux;
    acb_ptr roots;
    fmpz_t den;
    int res;
    slong k;
    slong lowprec = ACB_THETA_AGM_LOWPREC;

    fmpz_mat_init(D, 2*g, 2*g);
    fmpz_mat_init(R, 2*g, 2*g);
    acb_mat_init(aux, g, g);
    roots = _acb_vec_init(n);
    fmpz_init(den);

    res = accept_naive(tau, prec);
    if (res)
    {
	flint_printf("(all_const_sqr) Fall back to naive.\n");
	theta_naive(th2, tau, prec);
	goto exit;
    }
    
    res = is_balanced(aux, D, tau, prec);
    if (res)
    {
	flint_printf("(all_const_sqr) Fall back to newton.\n");
	acb_theta_newton_all_const_sqr(th2, tau, prec);
	goto exit;
    }

    flint_printf("Before real reduction:\n"); acb_mat_printd(aux, 10);
    
    /* Reduce real part */
    acb_siegel_reduce_real(R, aux, prec);
    acb_siegel_transform(aux, R, aux, prec);

    flint_printf("After real reduction:\n"); acb_mat_printd(aux, 10);
        
    /* Compute th2 recursively, act by R^-1 */
    acb_theta_all_const_sqr(th2, aux, prec);

    flint_printf("At exit of recursive call:\n");
    for (k = 0; k < n*n; k++)
    {
	acb_printd(&th2[k], 10); flint_printf("\n");
    }
    flint_printf("\n");
    
    fmpz_mat_inv(R, den, R);
    acb_theta_transform_all_sqr_proj(th2, th2, R, prec);
    acb_siegel_transform(aux, R, aux, prec);
    
    flint_printf("After reverse real reduction:\n");
    for (k = 0; k < n; k++)
    {
	acb_printd(&th2[k], 10); flint_printf("\n");
    }
    flint_printf("\n");
    
    /* Attempt duplication at D.aux: get roots at low precision */
    res = lowprec_roots(roots, aux, D, lowprec);
    
    if (!res)
    {
	/* Some roots are zero: abort with naive algorithm */
	theta_naive(th2, tau, prec);
	goto exit;
    }

    /* Duplicate */
    acb_theta_transform_sqr_proj(th2, th2, D, prec);
    
    flint_printf("Before duplication:\n");
    for (k = 0; k < n; k++)
    {
	acb_printd(&th2[k], 10); flint_printf("\n");
    }
    flint_printf("\n");
    
    flint_printf("Roots:\n");
    for (k = 0; k < n; k++)
    {
	acb_printd(&roots[k], 10); flint_printf("\n");
    }
    flint_printf("\n");
    
    for (k = 0; k < n; k++)
    {
	acb_theta_agm_sqrt_lowprec(&th2[k], &th2[k], &roots[k], prec); 
    }
    acb_theta_dupl_all_const(th2, th2, g, prec);
    
    flint_printf("After duplication:\n");
    for (k = 0; k < n*n; k++)
    {
	acb_printd(&th2[k], 10); flint_printf("\n");
    }
    flint_printf("\n");
    
    /* Act by inverse of D and x2 */
    fmpz_mat_inv(D, den, D);
    acb_theta_transform_all_sqr_proj(th2, th2, D, prec);
    _acb_vec_scalar_mul_2exp_si(th2, th2, n*n, 1);

    
    flint_printf("After final transformation:\n");
    for (k = 0; k < n*n; k++)
    {
	acb_printd(&th2[k], 10); flint_printf("\n");
    }
    flint_printf("\n");

    goto exit;

exit:
    {
	fmpz_mat_clear(D);
	fmpz_mat_clear(R);
	acb_mat_clear(aux);
	_acb_vec_clear(roots, n);
	fmpz_clear(den);
    }
}
