---
title: Dot products of integer vectors
functions: _fmpz_vec_dot_general _fmpz_vec_dot
---

`_fmpz_vec_dot` computes $\sum_i a_i b_i$ for two vectors of integers,
and `_fmpz_vec_dot_general` the more general $s \pm \sum_i a_i
b_{i}$ or $s \pm \sum_i a_i b_{n-1-i}$, which covers the inner loops of
matrix multiplication (`fmpz_mat`), classical polynomial multiplication
and division (`fmpz_poly`), and triangular solving. They give the same
result as a loop of `fmpz_addmul`, but run considerably faster when many
entries are small, which is the common case.

## Algorithm {#algorithm}

The function keeps three accumulators: a signed three-word integer for
the products of two inline entries, and two limb arrays for the sums of
the positive and the negative products that involve a large entry.

- A product of two inline entries has at most $2\bits - 4$ bits and is
  added to the three-word accumulator with a few word operations; the
  accumulator cannot overflow for any vector that fits in memory.
- A product of a large and an inline entry is added to the positive or
  negative limb array with one `mpn_addmul_1`.
- A product of two large entries is formed with `flint_mpn_mul` in a
  temporary array and added with `mpn_add`.

Keeping positive and negative sums apart avoids sign changes of the
accumulator in the loop. The arrays start on the stack with room for 64
limbs, so nothing is allocated unless the sums grow larger. At the end the
three sums and the initial value are combined and the result is
demoted. Zero entries are skipped. The two inputs may be any vectors,
and `res` may be the same object as `initial`, but not an entry of
either vector.

## Example {#example}

A dot product of small entries with factorials $(5i)!$, of which the
last ones are large, in normal and reversed order and with an initial
value:

```c example name=fmpz_vec_dot_variants
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_vec.h>

int main(void)
{
    slong i, n = 6;
    fmpz * a = _fmpz_vec_init(n);
    fmpz * b = _fmpz_vec_init(n);
    fmpz_t s, t;

    fmpz_init(s);
    fmpz_init(t);

    for (i = 0; i < n; i++)
    {
        fmpz_set_si(a + i, (i % 2) ? -(i + 1) : (i + 1));   /* 1 -2 3 -4 5 -6 */
        fmpz_fac_ui(b + i, 5 * i);                          /* (5i)! */
    }

    _fmpz_vec_dot(s, a, b, n);
    flint_printf("sum a_i b_i = "); fmpz_print(s); flint_printf("\n");

    /* the same with the second vector reversed, as in a convolution */
    _fmpz_vec_dot_general(s, NULL, 0, a, b, 1, n);
    flint_printf("sum a_i b_(n-1-i) = "); fmpz_print(s); flint_printf("\n");

    /* t - sum a_i b_i, accumulating into the initial value */
    fmpz_set_ui(t, 1000);
    _fmpz_vec_dot_general(t, t, 1, a, b, 0, n);
    flint_printf("1000 - sum a_i b_i = "); fmpz_print(t); flint_printf("\n");

    fmpz_clear(s);
    fmpz_clear(t);
    _fmpz_vec_clear(a, n);
    _fmpz_vec_clear(b, n);
    return 0;
}
```
```output
sum a_i b_i = -93067248095481105707385839
sum a_i b_(n-1-i) = 15511205177530892639309394
1000 - sum a_i b_i = 93067248095481105707386839
```
