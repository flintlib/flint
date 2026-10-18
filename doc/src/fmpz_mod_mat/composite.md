---
title: Prime and composite moduli
summary: Which `fmpz_mod_mat` functions need a prime modulus, what the others do when the modulus is composite, and which functions work for every modulus.
order: 1
---

Linear algebra over $\Z/n\Z$ is simplest when $n$ is prime and
$\Z/n\Z$ is a field: every nonzero pivot can be inverted, and Gaussian
elimination computes ranks, inverses, solutions and echelon forms. For
composite $n$, a nonzero entry may be a zero divisor, rank is no longer
well defined in the usual sense, and elimination needs other algorithms.
`fmpz_mod_mat` never tests whether $n$ is prime; the functions fall into
three groups.

## Functions for every modulus {#any}

These are correct for any $n \ge 1$:

- arithmetic: `fmpz_mod_mat_add`, `fmpz_mod_mat_mul`,
  `fmpz_mod_mat_pow_ui`, the matrix-vector products, `fmpz_mod_mat_trace`,
  and the conversions and random matrices;
- `fmpz_mod_mat_det` and `fmpz_mod_mat_charpoly`, which fall back to the
  division-free Berkowitz algorithm [@Ber1984] when elimination meets a
  non-invertible pivot;
- `fmpz_mod_mat_similarity`, which uses no division;
- `fmpz_mod_mat_howell_form` and `fmpz_mod_mat_strong_echelon_form`,
  which compute canonical forms of the row span over $\Z/n\Z$ with
  extended gcds instead of inverses [@How1986; @StoMul1998; @FieHof2014].
  The Howell form is the right replacement for the reduced row echelon
  form when $n$ is composite.

## Functions that need a prime modulus {#prime}

The other elimination functions are wrappers of generic algorithms for
fields over the `fmpz_mod` ring (`gr_mat_lu`, `gr_mat_rref_lu`,
`gr_mat_nonsingular_solve_lu`, ...). For a composite modulus they give
the correct result as long as every pivot they choose is invertible,
which is the common case for random matrices when the prime factors of
$n$ are large. When a pivot is a zero divisor, the generic algorithm
fails, and:

- `fmpz_mod_mat_lu`, `fmpz_mod_mat_rank`, `fmpz_mod_mat_rref`,
  `fmpz_mod_mat_nullspace`, `fmpz_mod_mat_solve_tril` and
  `fmpz_mod_mat_solve_triu` throw an exception;
- `fmpz_mod_mat_solve`, `fmpz_mod_mat_can_solve` and `fmpz_mod_mat_inv`
  return 0, as if the matrix were singular, even if it is invertible
  modulo $n$ (an inverse exists exactly when the determinant is a unit);
- `fmpz_mod_mat_minpoly` declares the ring a field and may return a
  wrong result.

Exceptions can be caught with the mechanism of `flint_set_throw`, but
for composite moduli it is better to use the functions of the first
group, or to factor $n$ and work modulo each prime power.

## Example {#example}

The matrix below has determinant 2, a unit modulo 15, but its first
pivot 3 is a zero divisor.

```c example name=fmpz_mod_mat_composite
#include <flint/fmpz.h>
#include <flint/fmpz_mod.h>
#include <flint/fmpz_mod_mat.h>
#include <flint/fmpz_mod_poly.h>

int main(void)
{
    fmpz_mod_ctx_t ctx;
    fmpz_mod_mat_t A, B, H;
    fmpz_mod_poly_t p;
    fmpz_t d;
    slong rank;

    fmpz_init(d);
    fmpz_mod_ctx_init_ui(ctx, 15);
    fmpz_mod_mat_init(A, 2, 2, ctx);
    fmpz_mod_mat_init(B, 2, 2, ctx);
    fmpz_mod_poly_init(p, ctx);

    /* A = [[3, 1], [1, 1]] */
    fmpz_set_ui(fmpz_mod_mat_entry(A, 0, 0), 3);
    fmpz_set_ui(fmpz_mod_mat_entry(A, 0, 1), 1);
    fmpz_set_ui(fmpz_mod_mat_entry(A, 1, 0), 1);
    fmpz_set_ui(fmpz_mod_mat_entry(A, 1, 1), 1);

    /* division-free algorithms work for any modulus */
    fmpz_mod_mat_det(d, A, ctx);
    printf("det: "); fmpz_print(d); printf("\n");
    fmpz_mod_mat_charpoly(p, A, ctx);
    printf("charpoly: "); fmpz_mod_poly_print_pretty(p, "x", ctx); printf("\n");

    /* elimination for fields gives up at the pivot 3 */
    printf("inv: %d\n", fmpz_mod_mat_inv(B, A, ctx));

    /* the Howell form works for any modulus: the rows (3, 6) and
       (5, 10) generate the same module as (1, 2) */
    fmpz_mod_mat_init(H, 3, 2, ctx);
    fmpz_set_ui(fmpz_mod_mat_entry(H, 0, 0), 3);
    fmpz_set_ui(fmpz_mod_mat_entry(H, 0, 1), 6);
    fmpz_set_ui(fmpz_mod_mat_entry(H, 1, 0), 5);
    fmpz_set_ui(fmpz_mod_mat_entry(H, 1, 1), 10);
    rank = fmpz_mod_mat_howell_form(H, ctx);
    printf("howell form, %ld nonzero rows:\n", rank);
    fmpz_mod_mat_print_pretty(H, ctx);
    printf("\n");

    fmpz_mod_mat_clear(A, ctx);
    fmpz_mod_mat_clear(B, ctx);
    fmpz_mod_mat_clear(H, ctx);
    fmpz_mod_poly_clear(p, ctx);
    fmpz_mod_ctx_clear(ctx);
    fmpz_clear(d);
    return 0;
}
```

```output
det: 2
charpoly: x^2+11*x+2
inv: 0
howell form, 1 nonzero rows:
[[1 2]
[0 0]
[0 0]
]
```

The characteristic polynomial is $x^2 - 4x + 2$, printed with residues
in $[0, 15)$.
