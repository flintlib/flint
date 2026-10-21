---
title: Dot products
functions: _gr_vec_dot _gr_vec_dot_strided _gr_vec_dot_si
---

The dot product $c \pm \sum_i a_i b_i$ is the central vector operation of
the generic interface. Classical polynomial multiplication computes each
coefficient of a product as a reversed dot product (`_gr_vec_dot_rev`),
classical and blocked matrix multiplication compute each entry as a dot
product of a row and a column (`_gr_vec_dot_strided` or a dot product
with a transposed operand), and triangular solving, LU decomposition,
power series inversion and Newton iteration are built from the same
primitive. A ring that implements only `GR_METHOD_VEC_DOT` and
`GR_METHOD_VEC_DOT_REV` efficiently therefore speeds up most generic
algorithms ([](algorithms/generic/method-tables#levels)).

## Arguments {#arguments}

All variants compute
$$
r = c + \sum_{i=0}^{n-1} a_i b_{\sigma(i)}
\quad\text{or}\quad
r = c - \sum_{i=0}^{n-1} a_i b_{\sigma(i)}
$$
where $n$ = `len`, $a$ = `vec1`, $b$ = `vec2`, the second form is used
when `subtract` is nonzero, and $c$ = `initial` may be `NULL`, which
means $c = 0$. The index map is $\sigma(i) = i$ for `_gr_vec_dot`,
$\sigma(i) = n - 1 - i$ for `_gr_vec_dot_rev`, and for
`_gr_vec_dot_strided` the entries used are $a_{i s_1}$ and $b_{i s_2}$
for strides $s_1$, $s_2$ counted in elements, which may be negative.
The variants `_gr_vec_dot_si`, `_gr_vec_dot_ui` and `_gr_vec_dot_fmpz`
take $b$ as a C array of integers, which is how integer matrices and
polynomials act on vectors over other rings without converting the
integers first. For $n = 0$ the result is $c$ (or zero).

The output `res` may alias `initial` but not an entry of `vec1` or
`vec2`. Implementations may evaluate the sum in any order and with
delayed reductions, so in floating-point rings the rounding differs from
that of a loop of `gr_addmul` calls; ball arithmetic rings return a
correct enclosure in any case.

## Implementations {#implementations}

The default `gr_generic_vec_dot` multiplies and adds through the method
table. The strided default, for rings that override `VEC_DOT` but not
`VEC_DOT_STRIDED`, gathers the strided entries into temporary arrays by
shallow copies and calls `_gr_vec_dot`. Specialised versions include:

- $\Z/n\Z$ (`nmod`, `nmod8`, `nmod32`, `nmod_redc`): accumulation of
  products in two or three words with a single reduction at the end, as
  in `_nmod_vec_dot` ([](algorithms/linear-algebra/dot-products)).
- $\Z/n\Z$ for 2 to 16 limbs (`mpn_mod`): unreduced multi-limb
  accumulation with delayed reduction
  ([](algorithms/arithmetic/multiword-delayed-reduction)).
- $\Z$ (`fmpz`): `_fmpz_vec_dot_general`, which accumulates the
  products in limb arrays, separately for small terms and for large
  positive and negative terms, without intermediate `fmpz` values.
- Real and complex balls (`arb`, `acb`): `arb_dot` and `acb_dot`, which
  sum exactly in a fixed-point accumulator and add a single error bound
  [@Joh2019dot].
- Floating point (`arf`, `acf`, `nfloat`, `dfloat`) and finite fields
  (`fq`, `fq_nmod`, `fq_zech`): their own vector dot products.

## Example {#example}

Dot products modulo 7 with the four variants, for $a = (1, 2, 3, 4)$,
$b = (5, 6, 0, 2)$, $c = 3$ and the integer weights
$w = (1, -2, 3, -4)$:

```c example name=gr_vec_dot
#include <flint/flint.h>
#include <flint/gr.h>
#include <flint/gr_vec.h>

int main(void)
{
    gr_ctx_t ctx;
    gr_ptr a, b, c, res;
    slong n = 4, w[4] = {1, -2, 3, -4};
    int status = GR_SUCCESS;

    gr_ctx_init_nmod(ctx, 7);
    a = gr_heap_init_vec(n, ctx);
    b = gr_heap_init_vec(n, ctx);
    c = gr_heap_init(ctx);
    res = gr_heap_init(ctx);

    status |= _gr_vec_set_str(a, "[1, 2, 3, 4]", n, ctx);
    status |= _gr_vec_set_str(b, "[5, 6, 0, 2]", n, ctx);
    status |= gr_set_ui(c, 3, ctx);

    status |= _gr_vec_dot(res, NULL, 0, a, b, n, ctx);
    flint_printf("a.b           = %{gr}\n", res, ctx);
    status |= _gr_vec_dot(res, c, 1, a, b, n, ctx);
    flint_printf("3 - a.b       = %{gr}\n", res, ctx);
    status |= _gr_vec_dot_rev(res, NULL, 0, a, b, n, ctx);
    flint_printf("a.rev(b)      = %{gr}\n", res, ctx);
    status |= _gr_vec_dot_strided(res, NULL, 0, a, 2, b, 2, 2, ctx);
    flint_printf("a0b0 + a2b2   = %{gr}\n", res, ctx);
    status |= _gr_vec_dot_si(res, c, 0, a, w, n, ctx);
    flint_printf("3 + a.w       = %{gr}\n", res, ctx);

    gr_heap_clear_vec(a, n, ctx);
    gr_heap_clear_vec(b, n, ctx);
    gr_heap_clear(c, ctx);
    gr_heap_clear(res, ctx);
    gr_ctx_clear(ctx);
    flint_printf("status %d\n", status);
    return 0;
}
```
```output
a.b           = 4
3 - a.b       = 6
a.rev(b)      = 5
a0b0 + a2b2   = 5
3 + a.w       = 0
status 0
```

Here $a \cdot b = 25 \equiv 4$, the reversed product is
$1 \cdot 2 + 2 \cdot 0 + 3 \cdot 6 + 4 \cdot 5 = 40 \equiv 5$, and
$3 + a \cdot w = 3 - 10 \equiv 0$.
