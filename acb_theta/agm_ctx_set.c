
#include "acb_theta.h"

/* Candidates for symplectic matrices */

static void
fmpz_mat_Mi(fmpz_mat_t N, slong i)
{
    slong g = fmpz_mat_nrows(N)/2;

    fmpz_mat_one(N);
    fmpz_one(fmpz_mat_entry(N, i, i+g));
    fmpz_set_si(fmpz_mat_entry(N, i+g, i), -1);
    fmpz_zero(fmpz_mat_entry(N, i+g, i+g));
}

static void
fmpz_mat_Nij(fmpz_mat_t N, slong i, slong j)
{  
    slong g = fmpz_mat_nrows(N)/2;

    fmpz_mat_one(N);
    fmpz_one(fmpz_mat_entry(N, i, j+g));
    fmpz_one(fmpz_mat_entry(N, j, i+g));
    fmpz_set_si(fmpz_mat_entry(N, i+g, j), -1);
    fmpz_set_si(fmpz_mat_entry(N, j+g, i), -1);
    fmpz_zero(fmpz_mat_entry(N, i+g, i+g));
    fmpz_zero(fmpz_mat_entry(N, j+g, j+g));
}

static void
acb_theta_agm_ctx_candidates(fmpz_mat_struct* Ni, slong try, slong g)
{
    slong j, u, v, c;
    fmpz_mat_t J;
    flint_rand_t state;

    flint_randinit(state);
    fmpz_mat_init(J, 2*g, 2*g);

    fmpz_mat_J(J);
    
    /* Change state according to try */
    for (j = 0; j < try; j++) n_randint(state, 2);
  
    fmpz_mat_one(&Ni[0]);
    if (g == 1)
    {
        fmpz_mat_J(&Ni[1]);
    }
    else if (g == 2 && try == 0)
    {
        for (j = 1; j <= g; j++)
        {
            fmpz_mat_Mi(&Ni[j], j-1);
        }
        j = g+1;
        for (u = 0; u < g; u++)
        {
            for (v = u+1; v < g; v++)
            {
                fmpz_mat_Nij(&Ni[j], u, v);
                j++;
            }
        }
    }
    else
    {
        for (j = 1; j < (1<<g); j++)
        {
            /* (JM)^2 for random upper triangular M */
            fmpz_mat_one(&Ni[j]);
            for (u = 0; u < g; u++)
            {
                for (v = u; v < g; v++)
                {
                    c = n_randint(state, 3) - 1;
                    fmpz_set_si(fmpz_mat_entry(&Ni[j], u, v+g), c);
                    fmpz_set_si(fmpz_mat_entry(&Ni[j], v, u+g), c);
                }
            }
            fmpz_mat_mul(&Ni[j], J, &Ni[j]);
            fmpz_mat_mul(&Ni[j], &Ni[j], &Ni[j]);
        }
    }
    flint_randclear(state);
    fmpz_mat_clear(J);
}

/* Collect data for a single matrix */

static void
acb_theta_agm_ctx_set_roots(acb_theta_agm_ctx_t ctx, slong k, slong prec)
{
    acb_mat_t Ntau;
    acb_ptr Nz;
    acb_t scal;
    slong g = acb_theta_agm_ctx_g(ctx);
    slong n = 1<<g;
    slong nb_bad;
    slong i;
    slong lowprec = ACB_THETA_AGM_LOWPREC;
    int is_ext = acb_theta_agm_ctx_is_ext(ctx);
  
    acb_mat_init(Ntau, g, g);
    Nz = _acb_vec_init(2*g);
    acb_init(scal);

    /* Transform z, tau */
    if (is_ext)
    {
        acb_siegel_transform_ext(Nz, Ntau, acb_theta_agm_ctx_matrix(ctx, k),
                acb_theta_agm_ctx_z(ctx), acb_theta_agm_ctx_tau(ctx), prec);
    }
    else
    {
        acb_siegel_transform(Ntau, acb_theta_agm_ctx_matrix(ctx, k),
                acb_theta_agm_ctx_tau(ctx), prec);
    }

    /* Compute number of bad steps */
    if (is_ext)
    {
        nb_bad = acb_theta_agm_ext_nb_bad_steps(Nz, Ntau, prec);
    }
    else
    {
        nb_bad = acb_theta_agm_nb_bad_steps(Ntau, prec);
    }
    nb_bad = nb_bad + 1;
    acb_theta_agm_ctx_reset_steps(ctx, k, nb_bad);

    /* Set roots to low precision */
    for (i = 0; i < nb_bad; i++)
    {
        if (is_ext)
        {
            acb_theta_naive(acb_theta_agm_ctx_roots(ctx, k) + 2*n*i,
                    Nz, 2, Ntau, lowprec);
        }
        else
        {
            acb_theta_naive_const(acb_theta_agm_ctx_roots(ctx, k) + i*n,
                Ntau, lowprec);
        }
        acb_mat_scalar_mul_2exp_si(Ntau, Ntau, 1);
    }

    /* Correct rescaling of roots */
    if (is_ext)
    {
        acb_inv(scal, &acb_theta_agm_ctx_roots(ctx, k)[n], lowprec);
        _acb_vec_scalar_mul(acb_theta_agm_ctx_roots(ctx, k),
                acb_theta_agm_ctx_roots(ctx, k), 2*n*nb_bad, scal, lowprec);
        acb_inv(scal, &acb_theta_agm_ctx_roots(ctx, k)[0], lowprec);
        for (i = 0; i < nb_bad; i++)
        {
            _acb_vec_scalar_mul(acb_theta_agm_ctx_roots(ctx, k) + 2*n*i,
                    acb_theta_agm_ctx_roots(ctx, k) + 2*n*i, n, scal, lowprec);
            acb_sqrt(scal, scal, lowprec);
        }
    
        flint_printf("(ctx_roots) Sequence of th_0:\n");
        for (i = 0; i < nb_bad; i++)
        {
            acb_printd(&acb_theta_agm_ctx_roots(ctx, k)[2*i*n], 10);
            flint_printf("\n");
            acb_printd(&acb_theta_agm_ctx_roots(ctx, k)[2*i*n+n], 10);
            flint_printf("\n\n");
        }
    }
    else
    {
        acb_inv(scal, &acb_theta_agm_ctx_roots(ctx, k)[0], lowprec);
        _acb_vec_scalar_mul(acb_theta_agm_ctx_roots(ctx, k),
                acb_theta_agm_ctx_roots(ctx, k), n*nb_bad, scal, lowprec);
    }

    /* Set k2, ab, eps */
    acb_theta_agm_ctx_k2(ctx, k)
        = acb_theta_k2(acb_theta_agm_ctx_matrix(ctx, k));
    acb_theta_agm_ctx_ab(ctx, k)
        = acb_theta_transform_image_char(acb_theta_agm_ctx_eps(ctx, k), 0,
                acb_theta_agm_ctx_matrix(ctx, k));

    acb_mat_clear(Ntau);
    _acb_vec_clear(Nz, 2*g);
    acb_clear(scal);
}

/* Compute bounds for a single matrix once roots, th have been set */

static void
acb_theta_agm_ctx_set_bounds(acb_theta_agm_ctx_t ctx, slong k, slong prec)
{
    acb_ptr a;
    arb_t abs;
    arb_t m, M;
    arf_t err;
    fmpz_t exp;
    slong g = acb_theta_agm_ctx_g(ctx);
    slong nb_th;
    slong nb_bad;
    slong i, j;
    
    nb_th = 1<<g;
    if (acb_theta_agm_ctx_is_ext(ctx)) nb_th *= 2;

    a = _acb_vec_init(nb_th);
    arb_init(abs);
    arb_init(m);
    arb_init(M);
    arf_init(err);
    fmpz_init(exp);

    nb_bad = acb_theta_agm_ctx_nb_bad_steps(ctx, k);
    
    /* Compute mi, Mi */
    for (i = 0; i < nb_bad; i++)
    {
        arb_pos_inf(m);
        arb_zero(M);
        for (j = 0; j < nb_th; j++)
        {
            acb_abs(abs, &acb_theta_agm_ctx_roots(ctx, k)[i*nb_th + j], prec);
            arb_min(m, m, abs, prec);
            arb_max(M, M, abs, prec);
        }
        arb_sqr(m, m, prec);
        arb_get_lbound_arf(&acb_theta_agm_ctx_mi(ctx, k)[i], m, prec);
        arb_sqr(M, M, prec);
        arb_get_ubound_arf(&acb_theta_agm_ctx_Mi(ctx, k)[i], m, prec);        
    }

    /* Set a to beginning of good steps; deduce convergence rates */
    for (i = 0; i < nb_th; i++)
    {
        acb_sqr(&a[i], &acb_theta_agm_ctx_roots(ctx, k)[(nb_bad-1)*nb_th + i],
                prec);        
    }
    if (acb_theta_agm_ctx_is_ext(ctx))
    {
        acb_theta_agm_ext_conv_rate(acb_theta_agm_ctx_c_ext(ctx, k),
                acb_theta_agm_ctx_c(ctx, k), acb_theta_agm_ctx_e(ctx, k),
                a, g, prec);
    }
    else
    {
        acb_theta_agm_conv_rate(acb_theta_agm_ctx_c(ctx, k),
                acb_theta_agm_ctx_e(ctx, k), a, g, prec);
        arf_zero(acb_theta_agm_ctx_c_ext(ctx, k));
    }
    
    /* Deduce minf */
    arf_max(err, acb_theta_agm_ctx_c(ctx, k), acb_theta_agm_ctx_c_ext(ctx, k));
    arb_set_arf(m, err);
    arb_sub_si(m, m, 1, prec);
    arb_neg(m, m);
    arb_mul_arf(m, m, &acb_theta_agm_ctx_mi(ctx, k)[nb_bad-1], prec);
    arb_get_lbound_arf(acb_theta_agm_ctx_minf(ctx, k), m, prec);

    /* Compute rad: this assumes th has been set */    
    acb_theta_agm_radius(acb_theta_agm_ctx_rad(ctx, k),
            acb_theta_agm_ctx_mi(ctx, k),
            acb_theta_agm_ctx_Mi(ctx, k),
            acb_theta_agm_ctx_minf(ctx, k),
            acb_theta_agm_ctx_nb_bad_steps(ctx, k), prec);
    if (acb_theta_agm_ctx_is_ext(ctx))
    {
        acb_theta_dupl_transform_radius(acb_theta_agm_ctx_rad(ctx, k),
            acb_theta_agm_ctx_rad(ctx, k),
            acb_theta_agm_ctx_th(ctx),
            acb_theta_agm_ctx_matrix(ctx, k), prec);
    }
    else
    {        
        acb_theta_dupl_transform_radius_const(acb_theta_agm_ctx_rad(ctx, k),
            acb_theta_agm_ctx_rad(ctx, k),
            acb_theta_agm_ctx_th(ctx),
            acb_theta_agm_ctx_matrix(ctx, k), prec);
    }

    /* Compute min, max */
    if (acb_theta_agm_ctx_is_ext(ctx))
    {
        acb_theta_agm_ext_rel_err(err, acb_theta_agm_ctx_c_ext(ctx, k),
                acb_theta_agm_ctx_e(ctx, k), 1, prec);
        arb_one(abs);
        arb_add_error_arf(abs, err);
        fmpz_one(exp);
        fmpz_mul_2exp(exp, exp, acb_theta_agm_ctx_nb_bad_steps(ctx, k));
        
        arb_set_arf(m, &acb_theta_agm_ctx_Mi(ctx, k)[nb_bad-1]);
        arb_div_arf(m, m, acb_theta_agm_ctx_minf(ctx, k), prec);
        arb_mul(m, m, abs, prec);
        arb_pow_fmpz(m, m, exp, prec);
        arb_get_ubound_arf(acb_theta_agm_ctx_max(ctx, k), m, prec);

        arb_set_arf(m, acb_theta_agm_ctx_minf(ctx, k));
        arb_div_arf(m, m, &acb_theta_agm_ctx_Mi(ctx, k)[nb_bad-1], prec);
        arb_mul(m, m, abs, prec);
        arb_pow_fmpz(m, m, exp, prec);
        arb_get_lbound_arf(acb_theta_agm_ctx_min(ctx, k), m, prec);
    }
    else
    {        
        arf_set(acb_theta_agm_ctx_max(ctx, k),
                &acb_theta_agm_ctx_Mi(ctx, k)[nb_bad-1]);
        arf_set(acb_theta_agm_ctx_min(ctx, k),
                acb_theta_agm_ctx_minf(ctx, k));
    }

    flint_printf("(set_bounds) Current bounds with matrix\n");
    fmpz_mat_print_pretty(acb_theta_agm_ctx_matrix(ctx, k)); flint_printf("\n");
    flint_printf("%wd bad steps, min, max, rad:\n", nb_bad);
    arf_printd(acb_theta_agm_ctx_min(ctx, k), 10); flint_printf("\n");
    arf_printd(acb_theta_agm_ctx_max(ctx, k), 10); flint_printf("\n");
    arf_printd(acb_theta_agm_ctx_rad(ctx, k), 10); flint_printf("\n");
    
    _acb_vec_clear(a, nb_th);
    arb_clear(abs);
    arb_clear(m);
    arb_clear(M);
    arf_clear(err);
    fmpz_clear(exp);
}

/* Set B3 */

static int
acb_theta_agm_ctx_set_B3(acb_theta_agm_ctx_t ctx, slong prec)
{
    slong dim = acb_theta_agm_ctx_dim(ctx);
    arf_t B2;
    fmpz_t e;
    slong exp;
    arb_t eta;
    acb_ptr r;
    acb_mat_t fd, fdinv;
    arb_t norm, bound, test;
    slong lowprec = ACB_THETA_AGM_LOWPREC;
    int res;
        
    arf_init(B2);
    fmpz_init(e);
    arb_init(eta);
    r = _acb_vec_init(dim);
    acb_mat_init(fd, dim, dim);
    acb_mat_init(fdinv, dim, dim);
    arb_init(norm);
    arb_init(bound);
    arb_init(test);
        
    /* Evaluate finite difference */
    acb_theta_cauchy(B2, acb_theta_agm_ctx_rho(ctx),
            acb_theta_agm_ctx_M(ctx), 2, dim, lowprec);
    arf_frexp(B2, e, B2);    
    exp = fmpz_get_si(e);
    arf_mul_2exp_si(B2, B2, exp);
    
    arb_one(eta);
    arb_mul_2exp_si(eta, eta, FLINT_MIN(- exp - n_clog(dim, 2), - prec/2));
    acb_theta_newton_fd(r, fd, acb_theta_agm_ctx_th(ctx), eta, ctx, prec);

    flint_printf("Finite diff:\n");
    acb_mat_printd(fd, 10); flint_printf("\n");    
    res = acb_mat_inv(fdinv, fd, prec);

    if (!res) arb_pos_inf(norm);
    else acb_mat_ninf(norm, fdinv, lowprec);
    
    flint_printf("Inv, norm:\n");
    acb_mat_printd(fdinv, 10); flint_printf("\n");
    arb_printd(norm, 10); flint_printf("\n");
      
    /* Is ||FD^-1||*n*B2*eta less than 1? */        
    arb_mul_arf(bound, norm, B2, lowprec);
    arb_mul_si(bound, bound, dim, lowprec);
    arb_mul(bound, bound, eta, lowprec);
    arb_sub_si(test, bound, 1, lowprec);
    if (!arb_is_negative(test)) res = 0;

    arb_mul(bound, bound, norm, lowprec);
    arb_add(bound, bound, norm, lowprec);
    arb_get_ubound_arf(acb_theta_agm_ctx_B3(ctx), bound, lowprec);

    arf_clear(B2);
    fmpz_clear(e);
    arb_clear(eta);
    _acb_vec_clear(r, dim);
    acb_mat_clear(fd);
    acb_mat_clear(fdinv);
    arb_clear(norm);
    arb_clear(bound);
    arb_clear(test);
    return res;
}

/* Get logs */

static void
acb_theta_agm_ctx_set_logs(slong* log_th, slong* log_M, slong* log_rho,
        slong* log_B1, slong* log_B2, slong* log_B3,
        const acb_theta_agm_ctx_t ctx)
{
    arf_t c;
    fmpz_t e;
    slong dim = acb_theta_agm_ctx_dim(ctx);
    slong lowprec = ACB_THETA_AGM_LOWPREC;
    slong g = acb_theta_agm_ctx_g(ctx);
    slong nb_th = 1<<g;
    arb_t abs, m;
    slong k;

    if (acb_theta_agm_ctx_is_ext(ctx)) nb_th *= 2;
    arf_init(c);
    fmpz_init(e);
    arb_init(abs);
    arb_init(m);
  
    arf_frexp(c, e, acb_theta_agm_ctx_M(ctx));
    *log_M = fmpz_get_si(e);
    arf_frexp(c, e, acb_theta_agm_ctx_rho(ctx));
    *log_rho = fmpz_get_si(e) - 1;  
    arf_mul_si(c, acb_theta_agm_ctx_M(ctx), 2*(dim+1), lowprec, ARF_RND_CEIL);
    arf_div(c, c, acb_theta_agm_ctx_rho(ctx), lowprec, ARF_RND_CEIL);
    arf_frexp(c, e, c);
    *log_B1 = fmpz_get_si(e);
    arf_mul_si(c, acb_theta_agm_ctx_M(ctx), 2*(dim+1)*(dim+2), lowprec,
            ARF_RND_CEIL);
    arf_div(c, c, acb_theta_agm_ctx_rho(ctx), lowprec, ARF_RND_CEIL);
    arf_div(c, c, acb_theta_agm_ctx_rho(ctx), lowprec, ARF_RND_CEIL);
    arf_frexp(c, e, c);
    *log_B2 = fmpz_get_si(e);  
    arf_frexp(c, e, acb_theta_agm_ctx_B3(ctx));
    *log_B3 = fmpz_get_si(e);

    arb_zero(m);
    for (k = 0; k < nb_th; k++)
    {
        acb_abs(abs, &acb_theta_agm_ctx_th(ctx)[k], lowprec);
        arb_max(m, m, abs, lowprec);        
    }
    arb_get_ubound_arf(c, m, lowprec);
    arf_frexp(c, e, c);
    *log_th = fmpz_get_si(e);

    arf_clear(c);
    fmpz_clear(e);
    arb_clear(abs);
    arb_clear(m);
}

/* User function */

void
acb_theta_agm_ctx_set(acb_theta_agm_ctx_t ctx, slong prec)
{  
    acb_mat_t half;
    arf_t m;
    slong lowprec = ACB_THETA_AGM_LOWPREC;  
    slong n = acb_theta_agm_ctx_nb(ctx);
    slong g = acb_theta_agm_ctx_g(ctx);
    slong k;
    int res;
    slong try = -1;
    int is_ext = acb_theta_agm_ctx_is_ext(ctx);

    acb_mat_init(half, g, g);
    arf_init(m);

    /* Set theta values */
    acb_mat_scalar_mul_2exp_si(half, acb_theta_agm_ctx_tau(ctx), -1);
    if (is_ext)
    {
        acb_theta_naive_proj(acb_theta_agm_ctx_th(ctx),
                acb_theta_agm_ctx_z(ctx), 2, half, prec);
    }
    else
    {
        acb_theta_naive_const_proj(acb_theta_agm_ctx_th(ctx), half, prec);
    }
    
    /* Try different matrix setups, starting with try=0 */
    while (try < ACB_THETA_AGM_NB_MATRIX_SETUPS)
    {
        try++;        
        acb_theta_agm_ctx_candidates(acb_theta_agm_ctx_matrix(ctx, 0), try, g);
        arf_pos_inf(acb_theta_agm_ctx_rho(ctx));
        arf_zero(acb_theta_agm_ctx_M(ctx));
      
        for (k = 0; k < n; k++)
	{
            acb_theta_agm_ctx_set_roots(ctx, k, prec);
            acb_theta_agm_ctx_set_bounds(ctx, k, lowprec);
            
            arf_min(acb_theta_agm_ctx_rho(ctx),
                    acb_theta_agm_ctx_rho(ctx),
                    acb_theta_agm_ctx_rad(ctx, k));
            arf_div(m, acb_theta_agm_ctx_max(ctx, k),
                    acb_theta_agm_ctx_min(ctx, 0), lowprec, ARF_RND_CEIL);
            arf_max(acb_theta_agm_ctx_M(ctx),
                    acb_theta_agm_ctx_M(ctx), m);
	}
        
        flint_printf("(ctx_set) current M, rho:\n");
        arf_printd(acb_theta_agm_ctx_M(ctx), 10); flint_printf("\n");
        arf_printd(acb_theta_agm_ctx_rho(ctx), 10); flint_printf("\n");

        if (!arf_is_finite(acb_theta_agm_ctx_M(ctx))) continue;
        if (arf_cmp_si(acb_theta_agm_ctx_rho(ctx), 0) <= 0) continue;

        res = acb_theta_agm_ctx_set_B3(ctx, prec);
        if (res) break;

        flint_printf("Invalid B3:\n");
        arf_printd(acb_theta_agm_ctx_B3(ctx), 10); flint_printf("\n");
    }
    
    acb_theta_agm_ctx_set_logs(&acb_theta_agm_ctx_log_th(ctx),
            &acb_theta_agm_ctx_log_M(ctx), &acb_theta_agm_ctx_log_rho(ctx),
            &acb_theta_agm_ctx_log_B1(ctx), &acb_theta_agm_ctx_log_B2(ctx),
            &acb_theta_agm_ctx_log_B3(ctx), ctx);
    
    acb_mat_clear(half);
    arf_clear(m);
}
