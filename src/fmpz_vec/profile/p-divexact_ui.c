#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "profiler.h"

void
_fmpz_vec_scalar_divexact_ui_naive(fmpz * vec1, const fmpz * vec2,
                               slong len2, ulong c)
{
    slong i;
    for (i = 0; i < len2; i++)
        fmpz_divexact_ui(vec1 + i, vec2 + i, c);
}

int main(void)
{
    slong len, bits;

    flint_rand_t state;
    flint_rand_init(state);

    fmpz *A, *Ac, *B;
    ulong c;

    double t1, t2, FLINT_SET_BUT_UNUSED(__);

    flint_printf("   len    bits(A*c)  bits(c)      old       new    speedup\n");

    for (len = 1; len <= 10; len++)
    {
        for (bits = 5; bits <= 10000; bits *= 2)
        {
            A = _fmpz_vec_init(len);
            Ac = _fmpz_vec_init(len);
            B = _fmpz_vec_init(len);

            _fmpz_vec_randtest(A, state, len, bits);
            c = n_randtest_not_zero(state);
            _fmpz_vec_scalar_mul_ui(Ac, A, len, c);

            TIMEIT_START;
            _fmpz_vec_scalar_divexact_ui_naive(B, Ac, len, c);
            TIMEIT_STOP_VALUES(__, t1);
            TIMEIT_START;
            _fmpz_vec_scalar_divexact_ui(B, Ac, len, c);
            TIMEIT_STOP_VALUES(__, t2);

            flint_printf("%6wd    %6wd   %6wd    %8g   %8g   %.3f\n",
                    len, _fmpz_vec_max_bits(Ac, len), FLINT_BIT_COUNT(c),
                    t1, t2, t1 / t2);

            _fmpz_vec_clear(A, len);
            _fmpz_vec_clear(Ac, len);
            _fmpz_vec_clear(B, len);
        }
    }

    flint_rand_clear(state);
}

