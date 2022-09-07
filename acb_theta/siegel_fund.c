
#include "acb_theta.h"

static
void fmpz_mat_siegel_fund_g1(fmpz_mat_t mat, slong j)
{
    switch(j)
    {
    case 0:
	fmpz_mat_J(mat);
	break;
    default:
	flint_printf("fmpz_mat_siegel_fund: Error (invalid index: %d)\n", j);
	flint_abort();
    }
}

static
void fmpz_mat_siegel_fund_g2(fmpz_mat_t mat, slong j)
{
    slong g = 2;
    fmpz_mat_t a, b, c, d;

    fmpz_mat_init(a, g, g);
    fmpz_mat_init(b, g, g);
    fmpz_mat_init(c, g, g);
    fmpz_mat_init(d, g, g);

    if (j < 15)
    {
	fmpz_mat_zero(a);
	fmpz_mat_one(c);
	fmpz_mat_neg(b, c);
    }
    if (15 <= j && j < 17)
    {
	fmpz_mat_one(a);
	fmpz_mat_neg(b, a);
    }
    if (17 <= j && j < 19)
    {
	fmpz_mat_zero(b);
	fmpz_one(fmpz_mat_entry(c, 0, 0));
	fmpz_set_si(fmpz_mat_entry(c, 0, 1), -1);
	fmpz_set_si(fmpz_mat_entry(c, 1, 0), -1);
	fmpz_one(fmpz_mat_entry(c, 1, 1));
    }
  
    switch(j)
    {
    case 0:
	fmpz_mat_zero(d);
	break;
    case 1:
	fmpz_one(fmpz_mat_entry(d, 0, 0));
	break;
    case 2:
	fmpz_set_si(fmpz_mat_entry(d, 0, 0), -1);
	break;
    case 3:
	fmpz_one(fmpz_mat_entry(d, 1, 1));
	break;
    case 4:
	fmpz_set_si(fmpz_mat_entry(d, 1, 1), -1);
	break;
    case 5:
	fmpz_mat_one(d);
	break;
    case 6:
	fmpz_mat_one(d);
	fmpz_mat_neg(d, d);
	break;
    case 7:
	fmpz_set_si(fmpz_mat_entry(d, 0, 0), -1);
	fmpz_one(fmpz_mat_entry(d, 1, 1));
	break;
    case 8:
	fmpz_one(fmpz_mat_entry(d, 0, 0));
	fmpz_set_si(fmpz_mat_entry(d, 1, 1), -1);
	break;
    case 9:
	fmpz_one(fmpz_mat_entry(d, 0, 1));
	fmpz_one(fmpz_mat_entry(d, 1, 0));
	break;
    case 10:
	fmpz_set_si(fmpz_mat_entry(d, 0, 1), -1);
	fmpz_set_si(fmpz_mat_entry(d, 1, 0), -1);
	break;
    case 11:
	fmpz_one(fmpz_mat_entry(d, 0, 0));
	fmpz_one(fmpz_mat_entry(d, 0, 1));
	fmpz_one(fmpz_mat_entry(d, 1, 0));
	break;
    case 12:
	fmpz_set_si(fmpz_mat_entry(d, 0, 0), -1);
	fmpz_set_si(fmpz_mat_entry(d, 0, 1), -1);
	fmpz_set_si(fmpz_mat_entry(d, 1, 0), -1);
	break;
    case 13:
	fmpz_one(fmpz_mat_entry(d, 0, 1));
	fmpz_one(fmpz_mat_entry(d, 1, 0));
	fmpz_one(fmpz_mat_entry(d, 1, 1));
	break;
    case 14:
	fmpz_set_si(fmpz_mat_entry(d, 0, 1), -1);
	fmpz_set_si(fmpz_mat_entry(d, 1, 0), -1);
	fmpz_set_si(fmpz_mat_entry(d, 1, 1), -1);
	break;
    case 15:
	fmpz_one(fmpz_mat_entry(c, 0, 0));
	fmpz_one(fmpz_mat_entry(d, 1, 1));
	break;
    case 16:
	fmpz_one(fmpz_mat_entry(c, 1, 1));
	fmpz_one(fmpz_mat_entry(d, 0, 0));
	break;
    case 17:
	fmpz_mat_one(a);
	fmpz_mat_one(d);
	break;
    case 18:
	fmpz_mat_one(a);
	fmpz_mat_neg(a, a);
	fmpz_mat_set(d, a);
	break;
    default:
	flint_printf("fmpz_mat_siegel_fund: Error (invalid index: %d)\n", j);
	flint_abort();
    }

    fmpz_mat_set_abcd(mat, a, b, c, d);
  
    fmpz_mat_clear(a);
    fmpz_mat_clear(b);
    fmpz_mat_clear(c);
    fmpz_mat_clear(d);
}

void
fmpz_mat_siegel_fund(fmpz_mat_t mat, slong j)
{
    slong g = fmpz_mat_nrows(mat)/2;
    switch(g)
    {
    case 1:
	fmpz_mat_siegel_fund_g1(mat, j);
	break;
    case 2:
	fmpz_mat_siegel_fund_g2(mat, j);
	break;
    default:
	flint_printf("fmpz_mat_siegel_fund: Error (not implemented for g = %d)\n", g);
	flint_abort();
    }
}
