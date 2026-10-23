---
title: Elementary functions of balls
functions: mp_real_sin_cos_bits mp_real_tan_bits mp_real_sin_cos_pi_bits mp_real_sin_cos_pi_ui_div_ui mp_real_exp_bits mp_real_log1p_bits mp_real_exp_notab_log2
---

These functions evaluate exp, log, the trigonometric and hyperbolic
functions and their inverses on `mp_real_t` balls. They take a target
precision in bits rather than limbs and choose their working precision
themselves, so that the result has a relative accuracy of about
$2^{-\mathit{prec}}$ whenever the argument allows it.

## Common behaviour {#common}

- **Midpoint and radius.** The midpoint $m$ of the argument is
  evaluated exactly as given, with all its limbs, and the radius $\rho$
  is propagated at the end through a bound for the variation of the
  function over the ball (for exp, the factor $1 \pm (e^\rho - 1)$; for
  log, $\rho/(m - \rho)$). When the radius limits the accuracy, `prec`
  is first lowered to a few bits beyond what the radius allows, so that
  no time is spent on digits the radius would swamp.
- **Relative accuracy near zeros.** The results keep their relative
  accuracy where the function value is small: $\exp$ near 0 (through
  `mp_real_expm1_bits`), $\log$ near 1, $\sin$, $\tan$, $\operatorname{atan}$,
  $\sinh$, $\tanh$, $\operatorname{asin}$, $\operatorname{asinh}$ and
  $\operatorname{atanh}$ near 0, $\operatorname{acos}$ and
  $\operatorname{acosh}$ near 1. The argument reduction adds the leading
  zero bits of the reduced argument to the working precision. The
  functions of $\pi x$ reduce exactly and keep relative accuracy near
  every zero; $\sin x$ and $\cos x$ near nonzero multiples of $\pi$ show
  the cancellation of the reduction in their radius instead.
- **Domains.** The functions with a restricted domain return an `int`:
  1 on success, or 0 with the output set to zero when the ball is not
  certified to lie in the domain (a log of a ball containing 0, an
  asin of a ball reaching beyond $\pm 1$, a tangent of a ball that may
  contain a pole). The functions without a status throw where the
  result's exponent would leave the safe range: $\exp(y)$ for
  $|y| \ge 2^{\bits - 5}$.
- **Aliasing.** Outputs may alias the argument. Where a function has two
  outputs (`mp_real_sin_cos_bits`, `mp_real_sinh_cosh_bits`), either may
  be `NULL`.
- **Huge arguments.** The sine, cosine and tangent give $[0 \pm 1]$,
  resp. status 0, for $|x| \ge 2^{\max(65536,\, 4\,\mathit{prec})}$,
  where the reduction modulo $\pi/2$ would need more bits of $\pi$ than
  is reasonable.

The functions with explicit precisions in limbs, `mp_real_sin_cos_pi_ui_div_ui`,
`mp_real_tan_pi_ui_div_ui` and the table-free exponentials, follow the
same conventions with $B^{-n}$ in place of $2^{-\mathit{prec}}$.

## Algorithms {#algorithms}

Each function reduces the argument to a fixed-point number in $[0, 1)$
(by $\log 2$ for exp, by the exponent for log, by $\pi/2$ for the
trigonometric functions, exactly by $1/2$ for the functions of $\pi x$,
and to $|x| < 1$ for atan) and calls the kernel for the working size:
table-driven kernels up to three limbs, dense-table kernels up to 32
limbs, bitwise reduction up to 600 limbs, multi-prime (diophantine)
reduction up to 65536 limbs (Newton–Taylor inversion of it for log and
atan), and table-free bit-burst evaluation beyond. The hyperbolic and
inverse functions are compositions of exp, log and atan arranged to
avoid cancellation. The guide [](mp_real/fixed-point) lists the kernels,
and [](algorithms/floating-point/elementary-functions) describes them.

::: summary algorithms/floating-point/elementary-functions :::

## Example {#example}

```c example name=mp_real_elementary
#include <flint/flint.h>
#include <flint/mp_real.h>
#include <flint/arb.h>

static void show(const char * name, const mp_real_t x)
{
    arb_t a;
    arb_init(a);
    mp_real_get_arb(a, x);
    flint_printf("%s = ", name);
    arb_printn(a, 30, 0);
    flint_printf("\n");
    arb_clear(a);
}

int main(void)
{
    mp_real_t x, s, c, y;
    int ok;

    mp_real_init(x);
    mp_real_init(s);
    mp_real_init(c);
    mp_real_init(y);

    /* x = 3/2 exactly */
    mp_real_set_ui(x, 3);
    mp_real_mul_2exp_si(x, x, -1);

    mp_real_exp_bits(y, x, 100);
    show("exp(3/2)", y);
    mp_real_log_bits(y, x, 100);
    show("log(3/2)", y);
    mp_real_sin_cos_bits(s, c, x, 100);
    show("sin(3/2)", s);
    show("cos(3/2)", c);
    ok = mp_real_tan_bits(y, x, 100);
    flint_printf("tan status %d\n", ok);
    show("tan(3/2)", y);

    /* x = 1/2 + 2^-200: cos(pi x) keeps its relative accuracy */
    mp_real_set_ui(x, 1);
    mp_real_mul_2exp_si(x, x, -200);
    mp_real_set_ui(y, 1);
    mp_real_mul_2exp_si(y, y, -1);
    mp_real_add(x, x, y, 10);
    mp_real_sin_cos_pi_bits(s, c, x, 100);
    show("sin(pi x)", s);
    show("cos(pi x)", c);

    /* a ball containing 0 is outside the domain of log */
    mp_real_set_si(x, -1);
    mp_real_add_error_2exp_si(x, 1);
    ok = mp_real_log_bits(y, x, 100);
    flint_printf("log status %d\n", ok);

    mp_real_set_ui(x, 1);
    mp_real_mul_2exp_si(x, x, -1);
    mp_real_atanh_bits(y, x, 100);
    show("atanh(1/2)", y);

    mp_real_clear(x);
    mp_real_clear(s);
    mp_real_clear(c);
    mp_real_clear(y);
    flint_cleanup_master();
    return 0;
}
```

```output
exp(3/2) = [4.48168907033806482260205546012 +/- 7.25e-31]
log(3/2) = [0.405465108108164381978013115464 +/- 3.50e-31]
sin(3/2) = [0.997494986604054430941723371141 +/- 4.88e-31]
cos(3/2) = [0.0707372016677029100881898514343 +/- 3.13e-32]
tan status 1
tan(3/2) = [14.1014199471717193876460836520 +/- 1.23e-29]
sin(pi x) = [1.00000000000000000000000000000 +/- 4.64e-69]
cos(pi x) = [-1.95501790801056086750878985681e-60 +/- 9.48e-93]
log status 0
atanh(1/2) = [0.549306144334054845697622618461 +/- 2.63e-31]
```

The value $\cos(\pi(1/2 + 2^{-200})) = -\sin(\pi 2^{-200})$ is about
$-\pi \cdot 2^{-200}$ and is computed with a relative error of about
$2^{-106}$, although the argument itself has 200 bits.
