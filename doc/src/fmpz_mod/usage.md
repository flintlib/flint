---
title: Using fmpz_mod
summary: How to set up an `fmpz_mod_ctx_t`, which functions need canonical residues, how to invert and divide modulo a composite number, and how `fmpz_mod` relates to the generic rings and to `mpn_mod`.
order: 1
---

In `fmpz_mod`, an element of $\Z/n\Z$ is an ordinary `fmpz_t` holding a
residue in $[0, n)$, and the modulus travels separately in an
`fmpz_mod_ctx_t`. There is no special element type: residues are
initialised, cleared, copied, compared and printed with the `fmpz`
functions, and `fmpz_mod` only supplies the operations that depend on
$n$.

## Usage {#usage}

1. Create a context with `fmpz_mod_ctx_init` (or `fmpz_mod_ctx_init_ui`)
   and clear it with `fmpz_mod_ctx_clear` when done. Initialisation
   precomputes what the arithmetic for the size of $n$ needs, so a
   context should be reused rather than created per operation.
2. Bring integers into canonical form with `fmpz_mod_set_fmpz`,
   `fmpz_mod_set_ui` or `fmpz_mod_set_si`.
3. Compute with `fmpz_mod_add`, `fmpz_mod_sub`, `fmpz_mod_neg`,
   `fmpz_mod_mul`, `fmpz_mod_pow_ui` and so on. These require canonical
   operands and return canonical results. Variants such as
   `fmpz_mod_add_si` or `fmpz_mod_mul_fmpz` take one operand that may be
   any integer.
4. Compare with `fmpz_equal` and `fmpz_is_zero`, and read the modulus
   with `fmpz_mod_ctx_modulus`.

The canonical-form requirement is not checked (except in builds with
assertions enabled): a residue outside $[0, n)$ silently gives wrong
results. `fmpz_mod_is_canonical` tests it and `fmpz_mod_assert_canonical`
throws if it fails.

## Division modulo a composite number {#division}

`fmpz_mod_inv` throws an exception when its argument is not invertible.
When $n$ may be composite, test with `fmpz_mod_is_invertible` first, or
use `fmpz_mod_divides`, which solves $ax \equiv b \pmod n$ whenever a
solution exists, also for non-invertible $a$. `fmpz_mod_pow_fmpz` with a
negative exponent returns 0 instead of throwing.

The program below works modulo $n = 10^{40} + 1$, which is divisible by
$17$.

```c example name=fmpz_mod_usage
#include <flint/fmpz.h>
#include <flint/fmpz_mod.h>

int main(void)
{
    fmpz_mod_ctx_t ctx;
    fmpz_t n, a, b, c;

    fmpz_init(n);
    fmpz_init(a);
    fmpz_init(b);
    fmpz_init(c);

    /* n = 10^40 + 1, divisible by 10^8 + 1 = 17 * 5882353 */
    fmpz_set_str(n, "10000000000000000000000000000000000000001", 10);
    fmpz_mod_ctx_init(ctx, n);

    /* reduce arbitrary integers to canonical residues */
    fmpz_mod_set_si(a, -5, ctx);
    fmpz_set_str(b, "123456789012345678901234567890123456789012345", 10);
    fmpz_mod_set_fmpz(b, b, ctx);
    fmpz_print(a); printf("\n");
    fmpz_print(b); printf("\n");

    /* arithmetic on canonical residues */
    fmpz_mod_mul(c, a, b, ctx);
    fmpz_mod_add_si(c, c, 7, ctx);
    fmpz_print(c); printf("\n");

    /* 17 divides n, so it is not invertible; 10 is */
    fmpz_set_ui(a, 17);
    printf("17 invertible: %d\n", fmpz_mod_is_invertible(a, ctx));
    fmpz_set_ui(a, 10);
    fmpz_mod_inv(c, a, ctx);
    fmpz_print(c); printf("\n");

    /* 17 x = 34 has solutions although 17 is not invertible */
    fmpz_set_ui(a, 17);
    fmpz_set_ui(b, 34);
    printf("divides: %d, ", fmpz_mod_divides(c, b, a, ctx));
    fmpz_print(c); printf("\n");

    /* negative exponents need an invertible base */
    fmpz_set_si(b, -1);
    printf("pow: %d\n", fmpz_mod_pow_fmpz(c, a, b, ctx));

    fmpz_mod_ctx_clear(ctx);
    fmpz_clear(n);
    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(c);
    return 0;
}
```

```output
9999999999999999999999999999999999999996
6789012345678901234567890123456789000000
6054938271605493827160549382716055000011
17 invertible: 0
9000000000000000000000000000000000000001
divides: 1, 2
pow: 0
```

`fmpz_mod_divides` returns one solution; here the solutions of
$17x \equiv 34$ are $x \equiv 2$ modulo $n/17$.

## Speed {#speed}

The cost of an operation depends on the size of $n$ (see
`fmpz_mod_ctx_init`): one-word moduli use `nmod` arithmetic, two-limb
moduli an inline Barrett reduction, and larger moduli general `fmpz`
arithmetic and division. Since residues above $2^{62}$ are stored as
pointers to GMP integers, every operation on them goes through the
`fmpz` layer, which costs much more than the arithmetic itself for
moduli of a few limbs.

For polynomials and matrices modulo $n$ with 2 to 16 limbs, the
generic ring `mpn_mod` (residues as fixed-length limb arrays, delayed
reduction in all basecase algorithms) is typically faster; see
[](mpn_mod/usage) and [](algorithms/arithmetic/multiword-overview).

## Generic rings {#generic}

`gr_ctx_init_fmpz_mod` creates the generic ring $\Z/n\Z$ with `fmpz`
residues; its elements are interchangeable with `fmpz_mod` residues, and
its arithmetic calls this module.
