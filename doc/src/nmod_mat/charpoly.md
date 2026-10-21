---
title: Computing characteristic and minimal polynomials
functions: nmod_mat_charpoly nmod_mat_minpoly
---

The characteristic polynomial $\chi_M(x) = \det(xI - M)$ of a square
matrix $M$ is monic of degree $r$, its coefficients are, up to sign, the
sums of the principal minors, and its constant term is
$(-1)^r \det M$. The minimal polynomial $\mu_M$ is the monic polynomial of
least degree with $\mu_M(M) = 0$; it divides $\chi_M$ and has the same
irreducible factors. These functions are declared in `nmod_poly.h`, and
the output polynomial must be initialised with the modulus of $M$.

::: summary algorithms/linear-algebra/charpoly :::

## Characteristic polynomial {#charpoly}

`nmod_mat_charpoly` works for any modulus and chooses as follows.

- $r \ge 32$, $n \le 255$ and $n$ prime: Danilevsky's algorithm on a copy
  of $M$ with entries stored in bytes (through the `gr` context for
  `nmod8`).
- $r > 8$: `nmod_mat_charpoly_danilevsky`; if it fails, which can only
  happen for a composite modulus,
- otherwise: `nmod_mat_charpoly_berkowitz`.

`nmod_mat_charpoly_danilevsky` reduces a copy of $M$ to a block triangular
matrix with companion matrices on the diagonal by $O(r)$ similarity
transforms, each costing $O(r^2)$, so $O(r^3)$ operations in all
[@Dan1937]. Each step divides by a pivot. If the pivot is zero, the
algorithm looks for a nonzero entry to swap in, and if there is none it
splits off a companion block. If the pivot is a nonzero zero divisor, it
gives up and returns 0. For a prime modulus it always succeeds.

`nmod_mat_charpoly_berkowitz` uses no division at all, so it is correct
over any commutative ring, at the cost of $O(r^4)$ operations, computed as
dot products with a single reduction each [@Ber1984]. For $r \le 8$ it is
faster than Danilevsky's algorithm, and `nmod_mat_det` uses it for
$5 \le r \le 8$.

## Minimal polynomial {#minpoly}

`nmod_mat_minpoly` requires a prime modulus. It takes unit vectors
$v = e_j$, computes the Krylov sequence $v, Mv, M^2 v, \ldots$ until the
first linear dependency, which gives the minimal polynomial of $v$ with
respect to $M$, and replaces the current result by its least common
multiple with that polynomial. The next unit vector is chosen outside the
span of the Krylov vectors found so far, which is maintained in echelon
form with `nmod_mat_reduce_row`. The process stops when the Krylov
vectors span the whole space, and the least common multiple is then the
minimal polynomial of $M$. Each new vector costs a matrix-vector product
and the reduction of a row, $O(r^2)$ operations. For each generator the
number of vectors is one more than the degree of the polynomial found, so
the total is $r + 1$ when the first unit vector already generates the
whole space, as is often the case, and at most $r(r + 1)$. The cost is
therefore $O(r^3)$ in the first case and $O(r^4)$ in the worst case, and
each least common multiple costs $O(r^2)$ operations.

`nmod_mat_minpoly_with_gens` also marks with a 1 in the array `P` the
indices $j$ of the unit vectors used. These vectors generate the space as
a module over $\Z/p\Z[x]$ with $x$ acting as $M$.

## Moduli {#moduli}

For a composite modulus the characteristic polynomial is still well
defined and `nmod_mat_charpoly` computes it. The minimal polynomial is
not: the polynomials of least degree vanishing at $M$ need not be
unique, and `nmod_mat_minpoly` may raise a `FLINT_IMPINV` exception.

For $n = 1$, the characteristic polynomial is the zero polynomial.

## Example {#example}

A $4 \times 4$ matrix modulo 7 with characteristic polynomial
$(x - 2)^3 (x - 3)$ and minimal polynomial $(x - 2)^2 (x - 3)$, then a
$10 \times 10$ matrix modulo 12 on which Danilevsky's algorithm fails:

```c example name=nmod_mat_charpoly_minpoly
#include <flint/flint.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_poly.h>

int main(void)
{
    nmod_mat_t M, T;
    nmod_poly_t f, g;
    slong i;
    ulong gens[4] = { 0, 0, 0, 0 };

    /* M = diag(J, 3) where J is the 2 x 2 Jordan block for 2, and a
       further eigenvalue 2, modulo 7 */
    nmod_mat_init(M, 4, 4, 7);
    nmod_mat_entry(M, 0, 0) = 2; nmod_mat_entry(M, 0, 1) = 1;
    nmod_mat_entry(M, 1, 1) = 2;
    nmod_mat_entry(M, 2, 2) = 3;
    nmod_mat_entry(M, 3, 3) = 2;

    nmod_poly_init(f, 7);
    nmod_poly_init(g, 7);

    nmod_mat_charpoly(f, M);
    flint_printf("charpoly: "); nmod_poly_print_pretty(f, "x"); flint_printf("\n");

    nmod_mat_minpoly_with_gens(g, M, gens);
    flint_printf("minpoly:  "); nmod_poly_print_pretty(g, "x"); flint_printf("\n");
    flint_printf("generators e_j for j =");
    for (i = 0; i < 4; i++)
        if (gens[i])
            flint_printf(" %wd", i);
    flint_printf("\n");

    /* check that the minimal polynomial annihilates M */
    nmod_mat_init(T, 4, 4, 7);
    nmod_poly_evaluate_mat(T, g, M);
    flint_printf("minpoly(M) is zero: %d\n", nmod_mat_is_zero(T));

    /* a similarity transform does not change the characteristic polynomial */
    nmod_mat_similarity(M, 2, 5);
    nmod_mat_charpoly(g, M);
    flint_printf("charpoly unchanged: %d\n", nmod_poly_equal(f, g));

    nmod_mat_clear(M);
    nmod_mat_clear(T);
    nmod_poly_clear(f);
    nmod_poly_clear(g);

    /* modulo 12, Danilevsky's algorithm can fail; charpoly then falls
       back to Berkowitz's division-free algorithm */
    nmod_mat_init(M, 10, 10, 12);
    for (i = 0; i < 100; i++)
        nmod_mat_entry(M, i / 10, i % 10) = (i * i * i / 11 + (i / 10) * (i % 10)) % 12;
    nmod_poly_init(f, 12);
    nmod_poly_init(g, 12);
    flint_printf("danilevsky succeeded: %d\n", nmod_mat_charpoly_danilevsky(f, M));
    nmod_mat_charpoly(f, M);
    nmod_mat_charpoly_berkowitz(g, M);
    flint_printf("charpoly == berkowitz: %d, constant term %wu, det %wu\n",
        nmod_poly_equal(f, g), nmod_poly_get_coeff_ui(f, 0), nmod_mat_det(M));
    nmod_mat_clear(M);
    nmod_poly_clear(f);
    nmod_poly_clear(g);
    return 0;
}
```
```output
charpoly: x^4+5*x^3+2*x^2+5*x+3
minpoly:  x^3+2*x+2
generators e_j for j = 0 1 2 3
minpoly(M) is zero: 1
charpoly unchanged: 1
danilevsky succeeded: 0
charpoly == berkowitz: 1, constant term 5, det 5
```

Modulo 7, $(x - 2)^3 (x - 3) = x^4 + 5x^3 + 2x^2 + 5x + 3$. Each unit
vector is used as a generator, because each lies outside the span of the
Krylov vectors of the previous ones. Modulo 12, the constant term of the
characteristic polynomial of the $10 \times 10$ matrix equals its
determinant.
