
#include "acb_theta.h"

static void slong_vec_max(slong* r, slong* v1, slong* v2, slong d)
{
  slong k;
  for (k = 0; k < d; k++)
    {
      r[k] = FLINT_MAX(v1[k], v2[k]);
    }
}

void acb_theta_eld_fill(acb_theta_eld_t E, const arb_mat_t Y, const arb_t normsqr,
		  arb_srcptr offset, slong* last_coords,
		  ulong a, slong prec)
{
  slong k;
  slong min, mid, max, step;
  slong d = acb_theta_eld_dim(E);
  slong g = acb_theta_eld_ambient_dim(E);
  arf_t b;

  arf_init(b);

  /* Set input data */
  for (k = 0; k < g-d; k++)
    {
      E->last_coords[k] = last_coords[k];
    }
  _arb_vec_set(acb_theta_eld_offset(E), offset, d);
  arb_set(acb_theta_eld_normsqr(E), normsqr);

  /* Compute other data */
  arb_get_ubound_arf(b, normsqr, prec);
  arb_set_arf(acb_theta_eld_rad(E), b);
  if (!arb_is_positive(normsqr))
    {
      arb_zero(acb_theta_eld_rad(E));
    }
  else
    {
      arb_sqrt(acb_theta_eld_rad(E), acb_theta_eld_rad(E), prec);
      arb_div(acb_theta_eld_rad(E), acb_theta_eld_rad(E), arb_mat_entry(Y,d-1,d-1), prec);
    }
  
  arb_div(acb_theta_eld_ctr(E), &acb_theta_eld_offset(E)[d-1], arb_mat_entry(Y,d-1,d-1), prec);
  arb_neg(acb_theta_eld_ctr(E), acb_theta_eld_ctr(E));

  acb_theta_eld_interval(&min, &mid, &max, acb_theta_eld_ctr(E),
			 acb_theta_eld_rad(E), (a >> (g-d)) % 2, prec);
  step = 2;
  
  acb_theta_eld_min(E) = min;
  acb_theta_eld_mid(E) = mid;
  acb_theta_eld_max(E) = max;
  acb_theta_eld_step(E) = step;
  
  /* Get nb_pts, after induction on children if d>1 */
  if (min > max)
    {
      acb_theta_eld_nb_pts(E) = 0;
      for (k = 0; k < d; k++) acb_theta_eld_box(E, k) = 0;
    }
  else if (d == 1)
    {
      acb_theta_eld_nb_pts(E) = (max - min)/step + 1;
      acb_theta_eld_box(E, 0) = FLINT_MAX(max, -min);
    }
  else /* ((d > 1) && (min <= max)) */
    {  
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
      nr = (max - mid)/step + 1;
      nl = (mid - min)/step;
      acb_theta_eld_init_children(E, nr, nl);
      
      /* Set offset_mid, offset_diff */
      for (k = 0; k < d-1; k++)
	{
	  arb_set(&offset_diff[k], arb_mat_entry(Y, k, d-1));
	  arb_mul_si(&offset_mid[k], &offset_diff[k], mid, prec);
	  arb_mul_si(&offset_diff[k], &offset_diff[k], step, prec);
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
	  c = mid + k*step;
	  acb_theta_eld_next_normsqr(next_normsqr, normsqr, arb_mat_entry(Y,d-1,d-1),
				     acb_theta_eld_ctr(E), c, prec);
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
	  
	  c = mid - (k+1)*step;
	  acb_theta_eld_next_normsqr(next_normsqr, normsqr, arb_mat_entry(Y,d-1,d-1),
				     acb_theta_eld_ctr(E), c, prec);
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
  
  arf_clear(b);
}
