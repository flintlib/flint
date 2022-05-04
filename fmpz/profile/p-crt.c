#include "flint/fmpz.h"
#include "flint/profiler.h"

void
_fmpz_crt_combine(fmpz_t r1r2, fmpz_t m1m2, const fmpz_t r1, const fmpz_t m1, const fmpz_t r2, const fmpz_t m2)
{
    fmpz_invmod(m1m2, m1, m2);
    fmpz_mul(m1m2, m1m2, m1);
    fmpz_sub(r1r2, r2, r1);
    fmpz_mul(r1r2, r1r2, m1m2);
    fmpz_add(r1r2, r1r2, r1);
    fmpz_mul(m1m2, m1, m2);
    fmpz_mod(r1r2, r1r2, m1m2);
}

void
_fmpz_crt_combine_uiui(fmpz_t r1r2, fmpz_t m1m2, ulong r1, ulong m1, ulong r2, ulong m2)
{
    mp_limb_t M[2];

    /*
        M = m1 * m2
        c = invmod(m1, m2) * m1
        r1r2 = ((r2 - r1) * c + r1) mod M
        m1m2 = M
    */

    /* M = m1 * m2 */
    umul_ppmm(M[1], M[0], m1, m2);

    if (M[1] == 0)
    {
        mp_limb_t c, v;

        c = n_invmod(m1, m2) * m1;

        if (r2 >= r1)
            v = n_mulmod2(r2 - r1, c, M[0]);
        else
            v = n_mulmod2(n_negmod(r1 - r2, M[0]), c, M[0]);

        v = n_addmod(v, r1, M[0]);

        fmpz_set_ui(r1r2, v);
        fmpz_set_ui(m1m2, M[0]);
    }
    else
    {
        mp_limb_t c[2], t[4], q[3], r[3];

        umul_ppmm(c[1], c[0], n_invmod(m1, m2), m1);

        if (r2 >= r1)
        {
            t[2] = mpn_mul_1(t, c, 2, r2 - r1);

            mpn_add_1(t, t, 3, r1);
            mpn_tdiv_qr(q, r, 0, t, 3, M, 2);
        }
        else
        {
            sub_ddmmss(r[1], r[0], M[1], M[0], 0, r1 - r2);
            mpn_mul_n(t, c, r, 2);
            mpn_add_1(t, t, 4, r1);
            mpn_tdiv_qr(q, r, 0, t, 4, M, 2);
        }

        fmpz_set_uiui(r1r2, r[1], r[0]);
        fmpz_set_uiui(m1m2, M[1], M[0]);
    }
}

void
tree_crt(fmpz_t r, fmpz_t m, mp_srcptr residues, mp_srcptr primes, slong len)
{
    if (len == 0)
    {
        fmpz_zero(r);
        fmpz_one(m);
    }
    else if (len == 1)
    {
        fmpz_set_ui(r, residues[0]);
        fmpz_set_ui(m, primes[0]);
    }
    else if (len == 2)
    {
        _fmpz_crt_combine_uiui(r, m, residues[0], primes[0], residues[1], primes[1]);
    }
    else
    {
        fmpz_t r1, m1, r2, m2;

        fmpz_init(r1);
        fmpz_init(m1);
        fmpz_init(r2);
        fmpz_init(m2);

        tree_crt(r1, m1, residues, primes, len / 2);
        tree_crt(r2, m2, residues + len / 2, primes + len / 2, len - len / 2);
        _fmpz_crt_combine(r, m, r1, m1, r2, m2);

        fmpz_clear(r1);
        fmpz_clear(m1);
        fmpz_clear(r2);
        fmpz_clear(m2);
    }
}

void
fmpz_print1(const fmpz_t n)
{
//    printf("%lu", fmpz_fdiv_ui(n, 1000000000)); printf("\n");
}

void
benchmark(slong num_primes, slong prime_bits)
{
    flint_rand_t state;
    fmpz_comb_temp_t temp;
    fmpz_comb_t comb;
    mp_ptr primes, residues;
    fmpz_t res;
    slong k;

    flint_randinit(state);
    primes = flint_malloc(num_primes * sizeof(mp_limb_t));
    residues = flint_malloc(num_primes * sizeof(mp_limb_t));
    fmpz_init(res);

    primes[0] = n_nextprime(UWORD(1) << (prime_bits - 1), 0);
    for (k = 1; k < num_primes; k++)
        primes[k] = n_nextprime(primes[k-1], 0);

    for (k = 0; k < num_primes; k++)
        residues[k] = n_randint(state, primes[k]);

    printf("simple tree:   ");
    TIMEIT_START
    fmpz_t tmp;
    fmpz_init(tmp);
    tree_crt(res, tmp, residues, primes, num_primes);
    fmpz_clear(tmp);
    TIMEIT_STOP
    fmpz_print1(res);

    printf("multi CRT:     ");
    TIMEIT_START
    fmpz_comb_init(comb, primes, num_primes);
    fmpz_comb_temp_init(temp, comb);
    fmpz_multi_CRT_ui(res, residues, comb, temp, 0);
    fmpz_comb_clear(comb);
    fmpz_comb_temp_clear(temp);
    TIMEIT_STOP
    fmpz_print1(res);

    printf("multi precomp: ");
    fmpz_comb_init(comb, primes, num_primes);
    fmpz_comb_temp_init(temp, comb);
    TIMEIT_START
    fmpz_multi_CRT_ui(res, residues, comb, temp, 0);
    TIMEIT_STOP
    fmpz_comb_clear(comb);
    fmpz_comb_temp_clear(temp);
    fmpz_print1(res);

    flint_free(primes);
    flint_free(residues);
    fmpz_clear(res);
}

int main()
{
    slong len, bits;

    bits = 5;
    for (len = 1; len <= 4000; len = FLINT_MAX(len * 1.5, len + 1))
    {
        printf("bits = %ld, len = %ld\n", bits, len);
        benchmark(len, bits);
    }

    bits = 64;
    for (len = 1; len <= 4000; len = FLINT_MAX(len * 1.5, len + 1))
    {
        printf("bits = %ld, len = %ld\n", bits, len);
        benchmark(len, bits);
    }
}
