---
title: Fixed-point kernels
summary: The fixed-point layer of mp_real: number formats, the ulp error contract, and which elementary function kernel serves which precision.
order: 2
---

Below the ball level, `mp_real` works with fixed-point numbers in limb
arrays. These carry no exponent, sign or radius: the caller knows the
range of the value and receives a bound for the error in units of the
last place. The elementary function kernels, the Newton inverses and
roots and the floors of constants live at this level, and the ball
functions (`mp_real_exp_bits` and so on) are thin layers that reduce the
argument, call a kernel and turn its bound into a radius. Using the
kernels directly saves that overhead and gives control over the
algorithm, which is what `nfloat` and the tuning programs do.

## Formats {#formats}

A fixed-point number $(x, n)$ is an array of $n$ limbs representing
$\sum_{i<n} x_i B^{i-n} \in [0, 1)$; its ulp is $B^{-n}$. Outputs may
have one or two extra limbs above the fraction:

| Output | Limbs | Range | Functions |
|---|---|---|---|
| $(y, n)$ | $n$ fraction limbs | $[0, 1)$ | log1p, atan, atanh, $-\log$, constants below 1 |
| $(y, n + 1)$ | $n$ fraction limbs, units limb $y_n$ | $[0, B)$ | exp, sin, cos, tan, sinh, cosh, constants above 1 |
| $(q, n + 2)$ | $n$ fraction limbs, two integral limbs | $[0, B^2)$ | Newton inverses, quotients and roots |

The reduced-argument functions return $1 - \cos t$ rather than
$\cos t$, so that a value close to 1 does not lose its information to
the units limb.

## The error contract {#err}

A function with an ulp-accurate contract takes `ulong * err` right after
its outputs. It writes a bound $E$ with $|y - f(x)| \le E \cdot B^{-n}$,
where $x$ is the exact value of the input array; for two outputs the
bound holds for each. `err` may be `NULL`. Each function documents the
largest value it can write, and callers size their guard limbs from it:
a bound below $2^g$ ulps costs $g$ bits, so a caller wanting $p$
correct bits calls the kernel at $\lceil (p + g)/\bits \rceil$ limbs.
The bitwise kernels compute at the output precision and write an a
priori bound that grows with the reduction parameter (at most
$9r + 100$ for exp, with $r \le 768$); the kernels built on ball
arithmetic (reduced arguments, table-free, diophantine, Newton–Taylor)
write the bound computed at run time, usually a few ulps. Some bounds
are one-sided; this is stated where it holds.

The Newton inverses and roots (`_newton` suffix) are not ulp-accurate:
their bounds are relative to the result, like $4B^{-n}/a$ for $1/a$.

## Which kernel {#choice}

The ball functions choose the kernel by the number of limbs $n$ and by
the number of leading zero bits $z$ of the argument; the table gives the
choices on 64-bit machines.

| Limbs | exp | sin, cos, tan | log1p | atan |
|---|---|---|---|---|
| 1–3 | table-driven (`_mp_real_small_*`) | table-driven | table-driven | table-driven |
| 4–10 | dense tables (`_mp_real_exp_opt_4`, ...) | dense tables | dense tables | dense tables |
| 11–32 | run-time dense tables | run-time dense tables | run-time dense tables | run-time dense tables |
| 33–600 | `_mp_real_exp_bitwise_rs` | `_mp_real_sin_cos_bitwise_rs` | `_mp_real_log1p_bitwise_rs` | `_mp_real_atan_bitwise_rs` |
| 601–65536 | `_mp_real_exp_diophantine` | `_mp_real_sin_cos_diophantine` | `_mp_real_neglog_newton` | `_mp_real_atan_newton` |
| above | `_mp_real_exp_notab` | `_mp_real_sin_cos_notab` | `_mp_real_neglog_newton` | `_mp_real_atan_newton` |

For arguments with many leading zero bits, the series of
`_mp_real_exp_reduced` and `_mp_real_sin_cos_reduced` (or
`_mp_real_series_rs` and `_mp_real_series_tapered`) are used instead,
since their cost falls with $z$ while that of the table-based kernels
does not; the threshold depends on $n$, from about 20 bits at a few
limbs to thousands at 65536 limbs. The tables of the bitwise and
diophantine kernels are computed on first use and cached for the
process. A single evaluation at huge precision is cheaper without them:
`mp_real_exp_notab_squaring` and `mp_real_exp_agm` avoid all tables. See
[](algorithms/floating-point/elementary-functions) for the algorithms.

## Example {#example}

The program computes $\exp(1/10)$ at four limbs with three kernels. The
bitwise kernel's bound is its a priori $9r + 100$ for the depth $r = 34$
of the four-limb kernel; the diophantine kernel computes its bound at run
time. The output of a kernel converts to a ball with
`_mp_real_set_mpn_2exp` and `_mp_real_add_error_ulps_at`.

```c example name=mp_real_fixed_exp
#include <flint/flint.h>
#include <flint/mpn_extras.h>
#include <flint/mp_real.h>
#include <flint/arb.h>

#define N 4

int main(void)
{
    mp_real_t t;
    arb_t a;
    ulong x[N], y1[N + 1], y2[N + 1], y3[N + 1], d[N + 1], err, err1, err2, err3;

    mp_real_init(t);
    arb_init(a);

    /* x = 1/10 as a fixed-point number with N fraction limbs */
    mp_real_set_ui(t, 1);
    mp_real_div_ui(t, t, 10, N + 2);
    _mp_real_get_fixed(x, &err, t, N);
    flint_printf("x: err = %wu\n", err);

    /* exp(x) with N fraction limbs and a units limb, three ways */
    _mp_real_exp_bitwise_rs(y1, &err1, x, N, 0);
    _mp_real_exp_diophantine(y2, &err2, x, N);
    _mp_real_exp_notab(y3, &err3, x, N);
    flint_printf("bound: bitwise %wu, diophantine %wu, notab %wu\n",
        err1, err2, err3);

    /* the results differ by at most the sum of the bounds */
    if (mpn_cmp(y1, y2, N + 1) >= 0)
        mpn_sub_n(d, y1, y2, N + 1);
    else
        mpn_sub_n(d, y2, y1, N + 1);
    flint_printf("|bitwise - diophantine| <= %wu ulps: %d\n", err1 + err2,
        flint_mpn_zero_p(d + 1, N) && d[0] <= err1 + err2);

    /* the value: y1 B^-N as an exact ball, with err1 ulps added */
    _mp_real_set_mpn_2exp(t, y1, N + 1, -FLINT_BITS * N);
    _mp_real_add_error_ulps_at(t, (double) err1, -N);
    mp_real_get_arb(a, t);
    flint_printf("exp(x) = ");
    arb_printn(a, 40, 0);
    flint_printf("\n");

    mp_real_clear(t);
    arb_clear(a);
    flint_cleanup_master();
    return 0;
}
```

```output
x: err = 1
bound: bitwise 406, diophantine 3, notab 128
|bitwise - diophantine| <= 409 ulps: 1
exp(x) = [1.105170918075647624811707826490246668225 +/- 4.53e-40]
```

The input $1/10$ is itself rounded, with an error of one ulp, which a
caller adds to the result's error through the derivative ($e^{x} < 2$,
so two ulps).
