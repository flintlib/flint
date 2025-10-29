#include "fmpz.h"
#include "fmpz_mod.h"
#include "fmpz_mod_vec.h"

void _fmpz_mod_vec_rand(fmpz *A, flint_rand_t state, slong len, const fmpz_mod_ctx_t ctx)
{
    slong i;
    for (i = 0; i < len; i++)
        fmpz_randm(A + i, state, fmpz_mod_ctx_modulus(ctx));
}
