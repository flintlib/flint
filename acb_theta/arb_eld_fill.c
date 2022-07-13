
#include "acb_theta.h"

void arb_eld_fill(arb_eld_t E, const arb_mat_t Y, const arb_t normsqr,
		  arb_srcptr offset, slong* last_coords,
		  ulong a, slong prec)
{
  slong k;
  slong d = arb_eld_dim(E);
  slong g = arb_eld_ambient_dim(E);

  /* Set input data */
  for (k = 0; k < g-d; k++) E->last_coords[k] = last_coords[k];
  _arb_vec_set(arb_eld_offset(E), offset, d);
  arb_set(arb_eld_normsqr(E), normsqr);

  /* Compute other data */
  arb_sqrt(arb_eld_rad(E), normsqr, prec);
  arb_div(arb_eld_rad(E), arb_eld_rad(E), arb_mat_entry(Y,d-1,d-1), prec);
  
  arb_div(arb_eld_ctr(E), &arb_eld_offset(E)[d-1], arb_mat_entry(Y,d-1,d-1), prec);
  arb_neg(arb_eld_ctr(E), arb_eld_ctr(E));

  arb_eld_interval(&arb_eld_min(E), &arb_eld_mid(E), &arb_eld_max(E),
		   arb_eld_ctr(E), arb_eld_rad(E), (a >> (g-d)) % 2, prec);
  arb_eld_step(E) = 2;
  
  /* Get nb_pts, after induction on children if d>1 */
  if (arb_eld_min(E) > arb_eld_max(E))
    {
      arb_eld_nb_pts(E) = 0;
    }
  else if (d == 1)
    {
      arb_eld_nb_pts(E) = (arb_eld_max(E) - arb_eld_min(E))/arb_eld_step(E) + 1;
    }
  else /* ((d > 1) && (arb_eld_min(E) <= arb_eld_max(E))) */
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
      nr = (arb_eld_max(E) - arb_eld_mid(E))/arb_eld_step(E) + 1;
      nl = (arb_eld_mid(E) - arb_eld_min(E))/arb_eld_step(E);
      arb_eld_init_children(E, nr, nl);
      
      /* Set offset_mid, offset_diff */
      for (k = 0; k < d-1; k++)
	{
	  arb_set(&offset_diff[k], arb_mat_entry(Y, k, d-1));
	  arb_mul_si(&offset_mid[k], &offset_diff[k], arb_eld_mid(E), prec);
	  arb_mul_si(&offset_diff[k], &offset_diff[k], arb_eld_step(E), prec);
	}
      _arb_vec_add(offset_mid, offset_mid, offset, d-1, prec);
      for (k = 0; k < g-d; k++) next_coords[k+1] = last_coords[k];
      
      /* Set children recursively */
      arb_eld_nb_pts(E) = 0;
      _arb_vec_set(next_offset, offset_mid, d-1);
      for (k = 0; k < nr; k++)
	{
	  c = arb_eld_mid(E) + k*arb_eld_step(E);
	  arb_eld_next_normsqr(next_normsqr, normsqr, arb_mat_entry(Y,d-1,d-1),
			       arb_eld_ctr(E), c, prec);
	  next_coords[0] = c;
	  arb_eld_fill(arb_eld_rchild(E, k), Y, next_normsqr, next_offset,
		       next_coords, a, prec);
	  arb_eld_nb_pts(E) += arb_eld_nb_pts(arb_eld_rchild(E, k));
	  if (k < nr) _arb_vec_add(next_offset, next_offset, offset_diff, d-1, prec);
	}
      _arb_vec_set(next_offset, offset_mid, d-1);
      for (k = 0; k < nl; k++)
	{
	  _arb_vec_sub(next_offset, next_offset, offset_diff, d-1, prec);
	  c = arb_eld_mid(E) - (k+1)*arb_eld_step(E);
	  arb_eld_next_normsqr(next_normsqr, normsqr, arb_mat_entry(Y,d-1,d-1),
			       arb_eld_ctr(E), c, prec);
	  next_coords[0] = c;
	  arb_eld_fill(arb_eld_lchild(E, k), Y, next_normsqr, next_offset,
		       next_coords, a, prec);
	  arb_eld_nb_pts(E) += arb_eld_nb_pts(arb_eld_lchild(E, k));
	}

      arb_clear(next_normsqr);
      flint_free(next_coords);
      _arb_vec_clear(offset_diff, d-1);
      _arb_vec_clear(offset_mid, d-1);
      _arb_vec_clear(next_offset, d-1);
    }
}
