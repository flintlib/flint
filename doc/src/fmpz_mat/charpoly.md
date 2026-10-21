---
title: Characteristic and minimal polynomials
functions: fmpz_mat_charpoly fmpz_mat_minpoly
---

`fmpz_mat_charpoly` computes the characteristic polynomial
$\chi_A(x) = \det(xI - A)$ of a square integer matrix, and
`fmpz_mat_minpoly` its minimal polynomial, the monic polynomial $\mu_A$
of least degree with $\mu_A(A) = 0$. Both have integer coefficients, and
$\mu_A$ divides $\chi_A$. For all but tiny matrices both are computed
modulo word-size primes and reconstructed by the Chinese remainder
theorem; the results are always proven.

::: summary algorithms/linear-algebra/charpoly :::

## Characteristic polynomial {#charpoly}

`fmpz_mat_charpoly` chooses the method from the dimension $n$ and the bit
size $b$ of the largest entry:

| Condition | Method |
|---|---|
| $n \le 3$ | explicit formulas |
| $A = 0$ | $x^n$ |
| $n \le 8$, or $n < \lceil \log_2 (b + 1) \rceil$ | `fmpz_mat_charpoly_berkowitz` |
| otherwise | `fmpz_mat_charpoly_modular` |

`fmpz_mat_charpoly_berkowitz` is a wrapper of the generic
`_gr_mat_charpoly_berkowitz`, which uses $O(n^4)$ additions and
multiplications and no divisions
([](algorithms/linear-algebra/charpoly#berkowitz)); its intermediate
results are coefficients of characteristic polynomials of submatrices, so
their size stays controlled.

`fmpz_mat_charpoly_modular` computes the bound $C$ of
`fmpz_mat_charpoly_bound` on the absolute values of the coefficients,
chooses $\lceil (\log_2 C + 1)/63 \rceil$ primes above $2^{63}$, reduces
$A$ modulo all of them at once with `fmpz_mat_multi_mod_ui`, computes the
characteristic polynomials modulo the primes in parallel with
`nmod_mat_charpoly` ($O(n^\omega)$ or $O(n^3)$ operations each), and
reconstructs the coefficients in the symmetric range with
`_fmpz_vec_multi_CRT_ui`. Since the number of primes is fixed by the a
priori bound, there is no early termination. With
$\log_2 C \approx n(b + \log_2 n)$, about $n(b + \log_2 n)/63$
characteristic polynomials modulo primes are needed, each costing
$O(n^3)$ word operations, so the total is $O(n^4 (b + \log n))$ word
operations, plus the multimodular reduction and reconstruction.

The bound: the coefficient of $x^{n-k}$ is $\pm$ the sum of the
$\binom{n}{k}$ principal $k \times k$ minors of $A$, each of which is at
most the product of the $k$ largest row norms (and of the $k$ largest
column norms) by Hadamard's inequality.

## Minimal polynomial {#minpoly}

`fmpz_mat_minpoly` and `fmpz_mat_minpoly_modular` compute minimal
polynomials modulo primes $p$ above $2^{63}$ with
`nmod_mat_minpoly_with_gens`, which also records which standard basis
vectors were used as generators of Krylov subspaces. For each prime:

- if the degree modulo $p$ is smaller than the degree seen before, $p$ is
  discarded: the minimal polynomial modulo $p$ divides the reduction of
  $\mu_A$, so a smaller degree identifies a bad prime;
- if it is larger, all previous primes were bad and are discarded;
- otherwise the coefficients are combined by the Chinese remainder theorem.

When the reconstruction is unchanged by a new prime, the candidate $f$ is
tested: $f(A) e_i = 0$ is checked exactly for every recorded generator
$e_i$, by $\deg f$ matrix-vector products. The Krylov spaces of these
generators span $\Q^n$ (they span $\mathbb{F}_p^n$), so the test proves
$\mu_A \mid f$; and $\deg f$, the degree modulo a prime, is at most
$\deg \mu_A$. Hence $f = \mu_A$. A bound on the coefficients of $\mu_A$
[@Dum2007], from a bound $\rho$ on the spectral radius obtained with
ovals of Cassini, caps the number of primes.

## Related functions {#related}

- `nmod_mat_charpoly` and `nmod_mat_minpoly` work modulo a word-size
  modulus.
- `fmpq_mat_charpoly` and `fmpq_mat_minpoly` scale a rational matrix to
  an integer matrix by the common denominator $q$ of its entries, call
  these functions, and rescale the coefficient of $x^{n-k}$ by $q^{-k}$.
- `fmpz_poly_factor` factors the result, for example to find the
  rational eigenvalues.

## Example {#example}

A $10 \times 10$ matrix with a repeated eigenvalue structure: the
characteristic polynomial has degree 10 while the minimal polynomial has
degree 4. The Berkowitz and modular algorithms agree:

```c example name=fmpz_mat_charpoly_minpoly
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_mat.h>

int main(void)
{
    fmpz_mat_t A;
    fmpz_poly_t c, c2, m;
    slong i, n = 10;

    fmpz_mat_init(A, n, n);
    fmpz_poly_init(c);
    fmpz_poly_init(c2);
    fmpz_poly_init(m);

    /* block diagonal: eigenvalue 2 (blocks of sizes 2, 2, 1),
       eigenvalue -1 (blocks of sizes 2, 1, 1, 1) */
    for (i = 0; i < 5; i++)
        fmpz_set_si(fmpz_mat_entry(A, i, i), 2);
    for (i = 5; i < 10; i++)
        fmpz_set_si(fmpz_mat_entry(A, i, i), -1);
    fmpz_one(fmpz_mat_entry(A, 0, 1));
    fmpz_one(fmpz_mat_entry(A, 2, 3));
    fmpz_one(fmpz_mat_entry(A, 5, 6));

    /* conjugate by a unimodular matrix to hide the structure */
    {
        fmpz_mat_t U, Ui, T;
        fmpz_t den;
        flint_rand_t state;

        flint_rand_init(state);
        fmpz_init(den);
        fmpz_mat_init(U, n, n);
        fmpz_mat_init(Ui, n, n);
        fmpz_mat_init(T, n, n);
        fmpz_mat_one(U);
        fmpz_mat_randops(U, state, 30);
        fmpz_mat_inv(Ui, den, U);   /* den = +-1 */
        fmpz_mat_mul(T, U, A);
        fmpz_mat_mul(A, T, Ui);
        fmpz_mat_scalar_divexact_fmpz(A, A, den);
        fmpz_mat_clear(U);
        fmpz_mat_clear(Ui);
        fmpz_mat_clear(T);
        fmpz_clear(den);
        flint_rand_clear(state);
    }

    fmpz_mat_charpoly(c, A);
    flint_printf("charpoly: %{fmpz_poly}\n", c);
    fmpz_mat_charpoly_berkowitz(c2, A);
    flint_printf("berkowitz agrees: %d\n", fmpz_poly_equal(c, c2));
    fmpz_mat_charpoly_modular(c2, A);
    flint_printf("modular agrees: %d\n", fmpz_poly_equal(c, c2));

    fmpz_mat_minpoly(m, A);
    flint_printf("minpoly: %{fmpz_poly}\n", m);

    fmpz_mat_clear(A);
    fmpz_poly_clear(c);
    fmpz_poly_clear(c2);
    fmpz_poly_clear(m);
    return 0;
}
```
```output
charpoly: x^10 - 5 * x^9 + 30 * x^7 - 15 * x^6 - 81 * x^5 + 30 * x^4 + 120 * x^3 - 80 * x - 32
berkowitz agrees: 1
modular agrees: 1
minpoly: x^4 - 2 * x^3 - 3 * x^2 + 4 * x + 4
```

The characteristic polynomial is $(x - 2)^5 (x + 1)^5$ and the minimal
polynomial $(x - 2)^2 (x + 1)^2$, the exponents being the sizes of the
largest Jordan blocks.
