
#include "acb_theta.h"

static void
slong_vec_max(slong* r, slong* v1, slong* v2, slong d)
{
    slong k;
    for (k = 0; k < d; k++)
    {
        r[k] = FLINT_MAX(v1[k], v2[k]);
    }
}

static void
acb_theta_eld_next_normsqr(arb_t next_normsqr, const arb_t normsqr,
        const arb_t gamma, const arb_t v, slong k, slong prec)
{
    arb_t x;
    arb_init(x);

    /* Set next_normsqr to normsqr - (v + gamma*k)^2 */
    arb_mul_si(x, gamma, k, prec);
    arb_add(x, x, v, prec);
    arb_sqr(x, x, prec);
    arb_sub(next_normsqr, normsqr, x, prec);

    arb_clear(x);
}

static void
acb_theta_eld_init_children(acb_theta_eld_t E, slong nr, slong nl)
{
    slong d = acb_theta_eld_dim(E);
    slong g = acb_theta_eld_ambient_dim(E);
    slong k;

    if (nr > 0)
    {
        E->rchildren = flint_malloc(nr * sizeof(struct acb_theta_eld_struct));
        acb_theta_eld_nr(E) = nr;
        for (k = 0; k < nr; k++) acb_theta_eld_init(acb_theta_eld_rchild(E, k), d-1, g);
    }
    if (nl > 0)
    {
        E->lchildren = flint_malloc(nl * sizeof(struct acb_theta_eld_struct));
        acb_theta_eld_nl(E) = nl;
        for (k = 0; k < nl; k++) acb_theta_eld_init(acb_theta_eld_lchild(E, k), d-1, g);
    }
}

static void
acb_theta_eld_init_interval(acb_theta_eld E, const arb_mat_t Y, const arb_t normsqr,
        arb_srcptr offset, slong* last_coords, ulong a, slong prec)
{    
    slong min, mid, max;
    slong d = acb_theta_eld_dim(E);
    slong g = acb_theta_eld_ambient_dim(E);
    slong k;
    arf_t b;
    arb_t ctr, rad;
    
    arf_init(b);
    arb_init(ctr);
    arb_init(rad);

    for (k = 0; k < g-d; k++)
    {
        E->last_coords[k] = last_coords[k];
    }
    
    arb_get_ubound_arf(b, normsqr, prec);
    arb_set_arf(rad, b);
    if (!arb_is_positive(normsqr))
    {
        arb_zero(rad);
    }
    else
    {
        arb_sqrt(rad, rad, prec);
        arb_div(rad, rad, arb_mat_entry(Y,d-1,d-1), prec);
    }
  
    arb_div(ctr, &offset(E)[d-1], arb_mat_entry(Y,d-1,d-1), prec);
    arb_neg(ctr, ctr);

    acb_theta_eld_interval(&min, &mid, &max, ctr, rad, (a >> (g-d)) % 2, prec);
  
    acb_theta_eld_min(E) = min;
    acb_theta_eld_mid(E) = mid;
    acb_theta_eld_max(E) = max;
    
    arf_clear(b);
    acb_clear(ctr);
    arb_clear(rad);
}

/* Main recursive function */

static void
acb_theta_eld_fill_recursive(acb_theta_eld_t E, const arb_mat_t Y, const arb_t normsqr,
        arb_srcptr offset, slong* last_coords, ulong a, slong prec);

void
acb_theta_eld_fill(acb_theta_eld_t E, const arb_mat_t Y, const arb_t normsqr,
        arb_srcptr offset, slong* last_coords, ulong a, slong prec)
{
    slong min, max;
    slong d = acb_theta_eld_dim(E);
    slong g = acb_theta_eld_ambient_dim(E);

    acb_theta_eld_init_interval(E, Y, normsqr, offset, last_coords, a, prec);  
    min = acb_theta_eld_min(E);
    max = acb_theta_eld_max(E);
  
    /* Induction only if d > 1 and min <= max */
    if (min > max)
    {
        acb_theta_eld_nb_pts(E) = 0;
        for (k = 0; k < d; k++) acb_theta_eld_box(E, k) = 0;
    }
    else if (d == 1)
    {
        acb_theta_eld_nb_pts(E) = (max - min)/2 + 1;
        acb_theta_eld_box(E, 0) = FLINT_MAX(max, -min);
    }
    else
    {
        acb_theta_eld_fill_recursive(E, Y, normsqr, offset, last_coords, a, prec);
    }
}

static void
acb_theta_eld_fill_recursive(acb_theta_eld_t E, const arb_mat_t Y, const arb_t normsqr,
        arb_srcptr offset, slong* last_coords, ulong a, slong prec)
{
    slong d = acb_theta_eld_dim(E);
    slong g = acb_theta_eld_ambient_dim(E);
    slong min = acb_theta_eld_min(E);
    slong mid = acb_theta_eld_mid(E);
    slong max = acb_theta_eld_max(E);
    slong k;
    
    arb_t next_normsqr;
    slong* next_coords;
    arb_ptr offset_diff;
    arb_ptr offset_mid;
    arb_ptr next_offset;
    slong c;
    slong nr, nl;
    
    arb_init(next_normsqr);
    next_coords = flint_malloc((g-d+1) * sizeof(slong));
    offset_diff = _arb_vec_init(d-1);
    offset_mid = _arb_vec_init(d-1);
    next_offset = _arb_vec_init(d-1);
    
    /* Initialize children */
    nr = (max - mid)/2 + 1;
    nl = (mid - min)/2;
    acb_theta_eld_init_children(E, nr, nl);
      
    /* Set offset_mid, offset_diff */
    for (k = 0; k < d-1; k++)
    {
        arb_set(&offset_diff[k], arb_mat_entry(Y, k, d-1));
        arb_mul_si(&offset_mid[k], &offset_diff[k], mid, prec);
        arb_mul_si(&offset_diff[k], &offset_diff[k], 2, prec);
    }
    _arb_vec_add(offset_mid, offset_mid, offset, d-1, prec);
    for (k = 0; k < g-d; k++) next_coords[k+1] = last_coords[k];
    
    /* Set children recursively */
    acb_theta_eld_nb_pts(E) = 0;
    acb_theta_eld_box(E, d-1) = FLINT_MAX(max, -min);
    for (k = 0; k < d-1; k++) acb_theta_eld_box(E, k) = 0;
    
    _arb_vec_set(next_offset, offset_mid, d-1);
    for (k = 0; k < nr; k++)
    {
        c = mid + k*2;
        acb_theta_eld_next_normsqr(next_normsqr, normsqr, arb_mat_entry(Y,d-1,d-1),
                &offset[d-1], c, prec);
        next_coords[0] = c;
        acb_theta_eld_fill(acb_theta_eld_rchild(E, k), Y, next_normsqr, next_offset,
                next_coords, a, prec);
        
        acb_theta_eld_nb_pts(E) += acb_theta_eld_nb_pts(acb_theta_eld_rchild(E, k));
        slong_vec_max(E->box, E->box, acb_theta_eld_rchild(E,k)->box, d-1);
        if (k < nr) _arb_vec_add(next_offset, next_offset, offset_diff, d-1, prec);
    }
    
    _arb_vec_set(next_offset, offset_mid, d-1);
    for (k = 0; k < nl; k++)
    {
        _arb_vec_sub(next_offset, next_offset, offset_diff, d-1, prec);
	
        c = mid - (k+1)*2;
        acb_theta_eld_next_normsqr(next_normsqr, normsqr, arb_mat_entry(Y,d-1,d-1),
                &offset[d-1], c, prec);
        next_coords[0] = c;
        acb_theta_eld_fill(acb_theta_eld_lchild(E, k), Y, next_normsqr, next_offset,
                next_coords, a, prec);
        
        acb_theta_eld_nb_pts(E) += acb_theta_eld_nb_pts(acb_theta_eld_lchild(E, k));
        slong_vec_max(E->box, E->box, acb_theta_eld_lchild(E,k)->box, d-1);
    }
    
    arb_clear(next_normsqr);
    flint_free(next_coords);
    _arb_vec_clear(offset_diff, d-1);
    _arb_vec_clear(offset_mid, d-1);
    _arb_vec_clear(next_offset, d-1);
}
