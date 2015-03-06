#include "fmpz_mat.h"
#include "fmpz.h"

void
fmpz_mat_content(fmpz_t ret, const fmpz_mat_t A)
{
    slong i, j, k;
    int cmp;
    
    k = 0;
    fmpz_set_si(ret, k);

    if (A->r != 0 && A->c != 0)
    {
	for (i = 0; i < A->r; i++)
	{
	    for (j = 0; j < A->c; j++)
	    {
		fmpz_gcd(ret, ret, fmpz_mat_entry(A, i, j));
		cmp = fmpz_is_one(ret);
		if (cmp == 1)
		    break;
	    }
	    if (cmp == 1)
		break;
	}
    }
    
}