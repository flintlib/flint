
#include "acb_theta.h"

ulong
acb_theta_transform_image_char(fmpz_t eps, ulong ab, const fmpz_mat_t mat)
{
  slong g = fmpz_mat_nrows(mat)/2;
  fmpz_mat_t a, b, c, d;
  fmpz_mat_t mat_tp;
  fmpz_mat_t block; /* CD^t or AB^t */
  fmpz_mat_t alphabeta;
  fmpz_mat_t alpha, beta; /* These are windows, not initialized or freed */
  fmpz_mat_t Cvec_1, Cvec_2, Lvec;
  fmpz_mat_t coef;
    
  ulong res = 0;
  slong i;

  fmpz_mat_init(a, g, g);
  fmpz_mat_init(b, g, g);
  fmpz_mat_init(c, g, g);
  fmpz_mat_init(d, g, g);
  fmpz_mat_init(mat_tp, 2*g, 2*g);
  fmpz_mat_init(block, g, g);
  fmpz_mat_init(alphabeta, 2*g, 1);
  fmpz_mat_init(Cvec_1, g, 1);
  fmpz_mat_init(Cvec_2, g, 1);
  fmpz_mat_init(Lvec, 1, g);
  fmpz_mat_init(coef, 1, 1);

  fmpz_mat_get_a(a, mat);
  fmpz_mat_get_b(b, mat);
  fmpz_mat_get_c(c, mat);
  fmpz_mat_get_d(d, mat);
  fmpz_mat_transpose(mat_tp, mat);

  /* Compute blocks and substract diagonals in alphabeta */
  fmpz_mat_transpose(block, d);
  fmpz_mat_mul(block, c, block);
  for (i = 0; i < g; i++)
    {
      fmpz_sub(fmpz_mat_entry(alphabeta, i, 0),
	       fmpz_mat_entry(alphabeta, i, 0), fmpz_mat_entry(block, i, i));
    }
  fmpz_mat_transpose(block, b);
  fmpz_mat_mul(block, a, block);
  for (i = 0; i < g; i++)
    {
      fmpz_sub(fmpz_mat_entry(alphabeta, g+i, 0),
	       fmpz_mat_entry(alphabeta, g+i, 0), fmpz_mat_entry(block, i, i));
    }

  /* Turn ab into a 2g x 1 fmpz matrix, and update alphabeta */
  for (i = 0; i < 2*g; i++)
    {
      /* Least significant bits first */
      fmpz_add_si(fmpz_mat_entry(alphabeta, 2*g-1-i, 0),
		  fmpz_mat_entry(alphabeta, 2*g-1-i, 0), ab & 1);
      ab = ab >> 1;		  
    }

  /* Perform matrix-vector multiplication */
  fmpz_mat_mul(alphabeta, mat_tp, alphabeta);


  /* Compute eps */
  fmpz_mat_window_init(alpha, alphabeta, 0, 0, g, 1);
  fmpz_mat_window_init(beta, alphabeta, g, 0, 2*g, 1);

  fmpz_zero(eps);
  
  fmpz_mat_mul(Cvec_1, c, beta);
  fmpz_mat_mul(Cvec_2, b, alpha);
  fmpz_mat_transpose(Lvec, Cvec_2);
  fmpz_mat_mul(coef, Lvec, Cvec_1);
  fmpz_addmul_ui(eps, fmpz_mat_entry(coef, 0, 0), 2);

  fmpz_mat_mul(Cvec_1, b, alpha);
  fmpz_mat_mul(Cvec_2, d, alpha);
  fmpz_mat_transpose(Lvec, Cvec_2);
  fmpz_mat_mul(coef, Lvec, Cvec_1);
  fmpz_sub(eps, eps, fmpz_mat_entry(coef, 0, 0));

  fmpz_mat_mul(Cvec_1, a, beta);
  fmpz_mat_mul(Cvec_2, c, beta);
  fmpz_mat_transpose(Lvec, Cvec_2);
  fmpz_mat_mul(coef, Lvec, Cvec_1);
  fmpz_sub(eps, eps, fmpz_mat_entry(coef, 0, 0));

  fmpz_mat_transpose(block, b);
  fmpz_mat_mul(block, a, block);
  for (i = 0; i < g; i++)
    {
      fmpz_set(fmpz_mat_entry(Lvec, 0, i), fmpz_mat_entry(block, i, i));
    }
  fmpz_mat_mul(Cvec_1, d, alpha);
  fmpz_mat_mul(Cvec_2, c, beta);
  fmpz_mat_sub(Cvec_1, Cvec_1, Cvec_2);
  fmpz_mat_mul(coef, Lvec, Cvec_1);
  fmpz_addmul_ui(eps, fmpz_mat_entry(coef, 0, 0), 2);

  fmpz_mod_ui(eps, eps, 8); /* Formula involves zeta_8^eps */

  fmpz_mat_window_clear(alpha);
  fmpz_mat_window_clear(beta);
  
  /* Reduce alphabeta mod 2 & convert to ulong */
  for (i = 0; i < 2*g; i++)
    {
      res = res << 1;
      res += fmpz_tstbit(fmpz_mat_entry(alphabeta, i, 0), 0);
    }

  fmpz_mat_clear(mat_tp);
  fmpz_mat_clear(block);
  fmpz_mat_clear(alphabeta);
  fmpz_mat_clear(Cvec_1);
  fmpz_mat_clear(Cvec_2);
  fmpz_mat_clear(Lvec);
  fmpz_mat_clear(coef);

  return res;
}
