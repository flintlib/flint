---
title: Fixed-point transforms
functions: gr_dft_ctx_init_nfixed GR_DFT_NFIXED_MAX_NLIMBS gr_dft_nfixed_set_arf gr_dft_precomp_nfixed_bound
---

The fixed-point contexts let a transform run on multi-limb fixed-point
numbers, which are much cheaper than floating-point or ball arithmetic at
the same precision, while `gr_dft_precomp_nfixed_bound` gives rigorous
bounds for the magnitude of every intermediate value and for the error of
the output. The functions of [](gr_dft/acb) are built on them. They are
useful directly for numerical code that controls the scaling of its data
itself, for example to perform many transforms of the same shape with an
error bound computed once.

## Representation {#representation}

An element of the real context with $k$ limbs is a sign limb, 0 for
nonnegative and 1 for negative values, followed by $k$ limbs holding the
absolute value as a fraction in $[0, 1)$, least significant limb first.
The unit in the last place is $\mathrm{ulp} = 2^{-k \bits}$. An element
of the complex context is a real part followed by an imaginary part, so
it occupies $2(k + 1)$ limbs and has the layout that complex mode
requires. Zero may have either sign.

| Operation | Result |
|---|---|
| addition, subtraction, negation | exact; `GR_UNABLE` and a wrapped value if the result has absolute value $\ge 1$ |
| real multiplication | truncated toward zero, error below 1 ulp |
| complex multiplication | error per part at most 2 ulp for $k \le 2$, 4 ulp for $3 \le k < 22$ (four real products), 6 ulp for $k \ge 22$ (three real products) |
| `gr_div_ui` | truncated toward zero, error below 1 ulp per part |
| `gr_one`, `gr_neg_one`, `gr_set_si` | $\pm 1$ become $\pm(1 - \mathrm{ulp})$; other integers except 0 give `GR_DOMAIN` |

Overflow is detected without branches and reported through the status,
which the transforms accumulate. A transform that returns `GR_SUCCESS`
has therefore had no intermediate value outside $(-1, 1)$, whether or not
the input was scaled according to the bound. The contexts implement only
the operations needed by the transforms and the basic `gr` interface
(including `gr_equal`, `gr_is_zero`, printing and random generation);
`gr_is_one` and `gr_is_neg_one` always return `T_FALSE`.

## Plans {#plans}

A plan over the fixed-point contexts should be made with
`gr_dft_precomp_init_karatsuba(P, n, alg, flags, R, C)`, where `R` and `C`
are the real and complex contexts with the same number of limbs. This
gives complex mode, in which multiplications by $\pm 1$ and $\pm i$ are
exact; the three-multiplication tables of complex mode are never used in
fixed point, because their entries $d \pm c$ can exceed 1. Any algorithm
and length can be used. The root table is computed from $e^{-2\pi i/n}$
evaluated with `arb` at 64 extra bits, by products
$w^j = w^{\lfloor j/2 \rfloor} w^{\lceil j/2 \rceil}$ up to an eighth of
the circle and by symmetry beyond, which gives entries with errors of
$O(j)$ ulp. The error bound of the table is stored in the plan. The
inverse transform works too, with a truncating division by $n$.

## Scaling the input {#scaling}

The caller must scale the input so that no value reaches 1 in absolute
value. `gr_dft_precomp_nfixed_bound(&peak, &err, m, d, P)` gives, for
inputs of complex modulus at most $m$ and errors of at most $d$ ulp per
part, a bound `peak` for the modulus of every value in the transform and
a bound `err` for the output error per part, in ulp. Both are linear in
$m$ and $d$, so the procedure is:

1. compute `peak` for $m = 1$ and $d = 0$;
2. scale the input by a power of two $2^{-s}$ with
   $2^{-s} \cdot \mathit{peak} < 1$ and convert it with
   `gr_dft_nfixed_set_arf`, which adds an error below 1 ulp per part;
3. compute `err` for $m$ equal to the largest input modulus and $d = 1$;
4. after the transform, check the status, and convert the outputs back
   with `gr_dft_nfixed_get_arf` and multiplication by $2^s$; the error of
   each part is at most $\mathit{err} \cdot 2^{s - k \bits}$.

`peak` is at least about $n$, because the outputs can be $n$ times larger
than the inputs. The bounds are pessimistic, typically by one or two
orders of magnitude. They describe the forward transform and the inverse
transform without the division by $n$ (as computed by
`_gr_dft_precomp_raw`); the truncating division then adds less than 1 ulp
and divides the error by $n$. [](algorithms/fft/fixed-point-dft)
explains how the bounds are derived.

## Example {#example}

A transform of length 64 with one limb: the bound for unit inputs is
$\mathit{peak} \approx 64$, so the input is scaled by $2^{-7}$. The
observed error is far below the bound.

```c example name=gr_dft_nfixed_bound
#include <math.h>
#include <flint/flint.h>
#include <flint/arb.h>
#include <flint/acb.h>
#include <flint/gr.h>
#include <flint/gr_dft.h>

int main(void)
{
    gr_ctx_t R, C;
    gr_dft_pre_t P;
    slong j, n = 64, nlimbs = 1, rsz;
    double peak1, err1, peak, err, maxerr = 0.0;
    int e, status;
    gr_ptr x, y;
    acb_ptr u, v;
    arf_t t;

    GR_MUST_SUCCEED(gr_dft_ctx_init_nfixed(R, nlimbs));
    GR_MUST_SUCCEED(gr_dft_ctx_init_nfixed_complex(C, nlimbs));
    GR_MUST_SUCCEED(gr_dft_precomp_init_karatsuba(P, n, GR_DFT_ALG_AUTO, 0, R, C));
    rsz = R->sizeof_elem;

    /* bound for unit inputs, then choose the input scale 2^-e */
    gr_dft_precomp_nfixed_bound(&peak1, &err1, 1.0, 0.0, P);
    frexp(peak1, &e);                  /* peak1 < 2^e */
    gr_dft_precomp_nfixed_bound(&peak, &err, ldexp(1.0, -e), 1.0, P);
    flint_printf("peak(1) = %.4f, scale 2^-%d, peak = %.4f, error bound %.1f ulp\n",
        peak1, e, peak, err);

    x = gr_heap_init_vec(n, C);
    y = gr_heap_init_vec(n, C);
    u = _acb_vec_init(n);
    v = _acb_vec_init(n);
    arf_init(t);

    /* inputs of modulus at most 2^-e: ((j mod 5) - 2)/4 + (j mod 3)/4 i */
    for (j = 0; j < n; j++)
    {
        acb_set_si_si(u + j, (j % 5) - 2, j % 3);
        acb_mul_2exp_si(u + j, u + j, -2 - e);
        GR_MUST_SUCCEED(gr_dft_nfixed_set_arf(GR_ENTRY(x, 2 * j, rsz),
                        arb_midref(acb_realref(u + j)), R));
        GR_MUST_SUCCEED(gr_dft_nfixed_set_arf(GR_ENTRY(x, 2 * j + 1, rsz),
                        arb_midref(acb_imagref(u + j)), R));
    }

    status = gr_dft_precomp(y, x, P, C);
    flint_printf("status: %s\n", status == GR_SUCCESS ? "GR_SUCCESS" : "overflow");

    /* compare with the exact transform, in units of 2^-64 */
    gr_dft_acb(v, u, n, 256);
    for (j = 0; j < 2 * n; j++)
    {
        arb_ptr exact = (j % 2 == 0) ? acb_realref(v + j / 2) : acb_imagref(v + j / 2);
        gr_dft_nfixed_get_arf(t, GR_ENTRY(y, j, rsz), R);
        arf_sub(t, t, arb_midref(exact), ARF_PREC_EXACT, ARF_RND_DOWN);
        arf_mul_2exp_si(t, t, FLINT_BITS * nlimbs);
        maxerr = FLINT_MAX(maxerr, fabs(arf_get_d(t, ARF_RND_NEAR)));
    }
    flint_printf("observed error %.1f ulp: %s the bound\n", maxerr,
        maxerr <= err ? "within" : "exceeds");

    arf_clear(t);
    _acb_vec_clear(u, n);
    _acb_vec_clear(v, n);
    gr_heap_clear_vec(x, n, C);
    gr_heap_clear_vec(y, n, C);
    gr_dft_precomp_clear(P);
    gr_ctx_clear(C);
    gr_ctx_clear(R);
    return 0;
}
```
```output
peak(1) = 64.0000, scale 2^-7, peak = 0.5000, error bound 2650.0 ulp
status: GR_SUCCESS
observed error 35.0 ulp: within the bound
```

The program uses `FLINT_BITS` = 64 for the unit; the output on a 32-bit
machine differs.
