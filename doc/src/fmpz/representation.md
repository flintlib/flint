---
title: The fmpz representation
summary: How an fmpz stores small values inline and large values in an mpz, the canonical form, how FLINT caches mpz structures, and how to write fast code with the internal functions and GMP.
order: 2
---

An `fmpz` is a single word. Most integers that occur in practice, such as
polynomial coefficients over $\Z$, matrix entries before fraction-free
elimination, or exponents of floating-point numbers, fit in a word, and
FLINT stores them without any allocation; larger integers are kept in a
GMP `mpz` structure that the word points to. Ordinary use of the
module needs none of this guide: see [](fmpz/getting-started).

## The encoding {#encoding}

Let $B = 2^{\bits}$. An `fmpz` $w$ is an `slong`, and its two most
significant bits decide how it is read:

| Top bits | Meaning | Range |
|---|---|---|
| `00` | the value $w$ itself | $0 \le w \le 2^{\bits - 2} - 1$ |
| `11` | the value $w$ itself | $-(2^{\bits - 2} - 1) \le w < 0$, and $w = -2^{\bits-2}$ is never used |
| `01` | a pointer: $4 (w - 2^{\bits - 2})$ is the address of an `mpz_struct` | any value of the `mpz` |
| `10` | unused | |

The macro `COEFF_IS_MPZ(w)` tests for the tag `01`, `COEFF_TO_PTR(w)`
shifts the word left by two bits to recover the address, and
`PTR_TO_COEFF(p)` shifts an address right and sets the tag. No
information is lost because every `mpz_struct` is at least 4-byte
aligned. `COEFF_MAX` $= 2^{\bits - 2} - 1$ and `COEFF_MIN` $=
-$`COEFF_MAX` bound the inline values: $\pm(2^{62}-1)$ on 64-bit machines.

The symmetric range of $\bits - 1$ bits has two consequences. Negation
of an inline value never overflows, and the sum or difference of two
inline values always fits in an `slong`, so additions are done in a
machine word and the result is checked against `COEFF_MAX` afterwards. Products of two inline values need a double word
(`smul_ppmm`), after which `fmpz_set_signed_uiui` stores the result
inline or in an `mpz`.

The `mpz` behind a large value is an ordinary GMP integer (`_mp_alloc`,
`_mp_size` with the sign, `_mp_d` the limbs). FLINT's own headers access
it as `zz_struct` through `FMPZ_TO_ZZ`, so that `fmpz.h` does not need
`gmp.h`.

## Canonical form {#canonical}

Every `fmpz` passed to or returned by a public function is *canonical*: a
value that fits inline is stored inline, never in an `mpz`. Functions rely
on this. `fmpz_is_zero`, `fmpz_is_one` and `fmpz_is_pm1` compare the word
with a constant, `fmpz_equal` concludes that two words with different
tags hold different values, and the fast paths of arithmetic functions
test `COEFF_IS_MPZ` only.

A function that computes into an `mpz` must therefore move a result that
turns out to be small back inline before returning; this is called
*demotion*, and the opposite, giving a small `fmpz` an `mpz`, is
*promotion*. The internal functions are:

- `_fmpz_promote(f)`: returns the `mpz` of $f$, attaching one if $f$ is
  small; the value is not preserved. Used when $f$ is about to be
  overwritten by a result known to be large.
- `_fmpz_promote_val(f)`: the same, but the `mpz` receives the value of a
  small $f$.
- `_fmpz_demote_val(f)`: if the `mpz` of $f$ holds a value that fits
  inline, stores it inline and releases the `mpz`. This restores the
  canonical form after a GMP or `mpn` call whose result may be small.
- `_fmpz_demote(f)`: releases the `mpz` of $f$, if any, and sets $f$ to 0.
  Used before writing an inline value directly to the word.

`_fmpz_is_canonical` checks the invariant and is used in the test code.

## Memory management {#memory}

Allocating and initialising an `mpz` for every large value would make
loops over moderately large integers dominated by `malloc`. FLINT
therefore keeps a cache of initialised `mpz` structures. Its behaviour
depends on how FLINT was configured.

**Default build.** Each thread has its own free list (in thread-local
storage). When it is empty, `_fmpz_new_mpz` allocates a block of 16
memory pages, initialises every `mpz_struct` that fits in them with room
for `MPZ_MIN_ALLOC` $= 2$ limbs, and puts them on the free list.
`_fmpz_clear_mpz` returns an `mpz` to the free list of the calling thread
without freeing its limbs, so the next value of the same size needs no
allocation at all; an `mpz` with more than 64 limbs is first shrunk back
to 2 limbs, so that the cache does not hold on to large amounts of
memory. Each page starts with a pointer to the header of its block, which
records the thread that allocated it. If an `mpz` is cleared by another
thread, or after its thread has cleaned up, it is freed with `mpz_clear`
and counted in its block header, and the block is freed when all its
`mpz`s have been. `flint_cleanup` (and `flint_cleanup_master`) release
the cached `mpz`s of the calling thread through `_fmpz_cleanup`.

**`--enable-reentrant`.** Every `mpz` is allocated with `flint_malloc`
and `mpz_init2` and freed on clear; there is no cache, and
`_fmpz_cleanup` does nothing.

**`--with-gc`.** For use with the Boehm garbage collector: a single
global cache protected by a mutex, which also keeps pointers to all
allocated `mpz`s so that the collector does not reclaim them.

In all builds, an `mpz` handed out by `_fmpz_new_mpz` has at least
`MPZ_MIN_ALLOC` limbs, and an `mpz` handed back to `_fmpz_clear_mpz` must
still have that many. GMP functions never shrink an allocation, with one
known exception: `mpz_remove` in GMP 6.3 may reduce it. Code that calls
`mpz_realloc` on an `mpz` owned by an `fmpz` must not go below
`MPZ_MIN_ALLOC`.

The memory used by one integer is reported by `fmpz_allocated_bytes` from
[](fmpz_extras): 0 for an inline value, and the size of the `mpz_struct`
plus its limbs otherwise.

```c example name=fmpz_representation
#include <gmp.h>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_extras.h>

static void info(const char * name, const fmpz_t x)
{
    flint_printf("%s: %s, %wd bytes allocated, canonical %d\n", name,
        COEFF_IS_MPZ(*x) ? "mpz" : "inline",
        fmpz_allocated_bytes(x), _fmpz_is_canonical(x));
}

int main(void)
{
    fmpz_t x;
    mpz_t z;

    fmpz_init_set_si(x, COEFF_MAX);
    info("COEFF_MAX", x);
    fmpz_add_ui(x, x, 1);
    info("COEFF_MAX + 1", x);
    fmpz_sub_ui(x, x, 1);
    info("COEFF_MAX again", x);
    fmpz_mul_2exp(x, x, 1000);
    info("COEFF_MAX * 2^1000", x);

    /* read-only access for a GMP function, without copying */
    flint_mpz_init_set_readonly(z, x);
    flint_printf("popcount = %wu\n", (ulong) mpz_popcount(z));
    flint_mpz_clear_readonly(z);

    fmpz_clear(x);
    return 0;
}
```
```output
COEFF_MAX: inline, 0 bytes allocated, canonical 1
COEFF_MAX + 1: mpz, 32 bytes allocated, canonical 1
COEFF_MAX again: inline, 0 bytes allocated, canonical 1
COEFF_MAX * 2^1000: mpz, 152 bytes allocated, canonical 1
popcount = 62
```

$2^{62}$ needs an `mpz` of 16 bytes with the minimum of two limbs; after
the subtraction the value is demoted again. The shifted value has 1062
bits, 17 limbs, which the cached `mpz` grew to hold.

## Writing functions on the representation {#underscore}

Most FLINT functions on `fmpz` have the same shape: a fast path for
inline operands that works on words, and a general path that works on the
`mpz` (or its limbs) and demotes the result. New low-level code can
follow the same pattern. Three rules keep it correct.

1. Read the word once (`slong v = *x`) and test it with `COEFF_IS_MPZ`.
   Only write an inline value directly to an output word if the output is
   known to be inline (otherwise its `mpz` leaks), or after
   `_fmpz_demote`.
2. Before writing a result that may be large through GMP or `mpn`
   functions, obtain the output `mpz` with `_fmpz_promote` (or
   `_fmpz_promote_val` if its value is needed), and finish with
   `_fmpz_demote_val` unless the result cannot be small.
3. Handle aliasing: once an output has been promoted, an input that is
   the same object may have changed. FLINT functions usually compare
   pointers and fall back to a temporary.

The functions with a leading underscore that do these steps are listed
under [](fmpz#internal). For values that are almost always small, such as
exponents, the inline functions of [](fmpz_extras) (`fmpz_add_inline`,
`fmpz_add_si_inline`, `_fmpz_sub_small`) put the fast path in the
caller's code, and the macro `FMPZ_GET_MPN_READONLY` gives a limb pointer
for any value.

The following program takes the next prime after $10^{30}$ with a GMP
function applied to the `mpz` of an `fmpz`, then counts the steps of the
Collatz iteration $x \mapsto x/2$ or $3x + 1$ from it, with an inline
path for words.

```c example name=fmpz_fast_path
#include <gmp.h>
#include <flint/flint.h>
#include <flint/fmpz.h>

/* x -> x/2 or 3x + 1 for x > 0, with an inline path for small x */
static void collatz_step(fmpz_t x)
{
    slong v = *x;

    if (!COEFF_IS_MPZ(v) && v < COEFF_MAX / 3)
        *x = (v % 2 == 0) ? v / 2 : 3 * v + 1;  /* stays below COEFF_MAX */
    else if (fmpz_is_even(x))
        fmpz_fdiv_q_2exp(x, x, 1);
    else
    {
        fmpz_mul_ui(x, x, 3);
        fmpz_add_ui(x, x, 1);
    }
}

int main(void)
{
    fmpz_t x;
    slong steps;
    mpz_ptr z;

    fmpz_init(x);
    fmpz_set_str(x, "1000000000000000000000000000000", 10);

    /* call a GMP function on the mpz of x directly */
    z = _fmpz_promote_val(x);
    mpz_nextprime(z, z);
    _fmpz_demote_val(x);
    flint_printf("p = "); fmpz_print(x); flint_printf("\n");

    for (steps = 0; !fmpz_is_one(x); steps++)
        collatz_step(x);
    flint_printf("%wd steps to reach 1\n", steps);

    fmpz_clear(x);
    return 0;
}
```
```output
p = 1000000000000000000000000000057
764 steps to reach 1
```

Writing `*x` directly in `collatz_step` is safe because the branch is
taken only when $x$ is inline. When $x$ becomes small again,
`fmpz_fdiv_q_2exp` demotes it, and from then on every step is a few word
operations.

## Interfacing with GMP {#gmp}

There are three ways to hand an integer to GMP and back, in increasing
order of speed and care required.

- **Copies.** `fmpz_get_mpz` and `fmpz_set_mpz` convert by copying the
  limbs. Always safe.
- **Read-only views.** `flint_mpz_init_set_readonly(z, f)` makes an
  `mpz_t` that shares the limbs of a large $f$ (or holds a small $f$ in
  its own one-limb buffer), for passing to GMP functions that only read
  their argument; release it with `flint_mpz_clear_readonly`. In the other
  direction, `fmpz_init_set_readonly(f, z)` makes an `fmpz_t` view of an
  `mpz_t`, released with `fmpz_clear_readonly`. Neither side may be
  modified while the view exists.
- **Direct access.** `_fmpz_promote_val` returns the `mpz` of an `fmpz`,
  which any GMP function may then read and write, followed by
  `_fmpz_demote_val`, as in the example above. The `mpz` remains owned by
  FLINT: it must not be cleared with `mpz_clear`, and its allocation must
  not shrink below `MPZ_MIN_ALLOC` limbs.

All three require `gmp.h` to be included before the FLINT headers, since
`fmpz.h` declares the functions with `mpz_t` arguments only in that case.
