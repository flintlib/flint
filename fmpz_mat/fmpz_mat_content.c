#include "fmpz_mat.h"
#include "fmpz.h"

void
fmpz_mat_content(fmpz_t ret, const fmpz_mat_t A)
{
    slong i, j, k;
    fmpz_t t;
    
    k = 0;
    
    if(A->r == 0 && A->c == 0)
	fmpz_set_si(ret, k);
    else
    {
	ret = fmpz_mat_entry(A, 0, 0);
    
	for(i = 0; i < A->r; i++)
	{
	    for(j = 0; j < A->c; j++)
	    {
		t = ret;
		fmpz_gcd(ret, t, fmpz_mat_entry(A, i, j)); 
	    }
	}
    }
    
}