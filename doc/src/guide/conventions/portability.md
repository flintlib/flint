---
title: Portability
summary: Writing code with FLINT that works on 32-bit and 64-bit platforms and on Windows: word types, constants, printing and the platform-dependent parts of the library.
order: 4
---

FLINT runs on 64-bit and 32-bit platforms, on Unix-like systems and on
Windows. Code using FLINT stays portable if it follows a few rules about
word sizes.

## Word types and constants {#types}

The C types `long` and `unsigned long` have 64 bits on 64-bit Unix
systems but only 32 bits on 64-bit Windows. FLINT therefore uses its own
types `ulong` and `slong`, which always have the size of a GMP limb,
`FLINT_BITS` bits: 64 on all 64-bit platforms, including Windows, and 32
on 32-bit platforms. Use these types for all values passed to FLINT
functions taking `ulong` or `slong`, and do not assume that they are
`unsigned long` and `long`.

Integer constants of word size are written with the macros `UWORD` and
`WORD`, as in `UWORD(1) << 40` or `WORD(-123456789012)`, instead of the
suffixes `UL` and `L`, which give 32-bit constants on 64-bit Windows. The
extreme values are `UWORD_MAX`, `WORD_MIN` and `WORD_MAX`. Code that
depends on the word size can test `FLINT_BITS` (or `FLINT64`) with the
preprocessor, and the same holds for constants that only fit in 64 bits.

## Printing and reading words {#printing}

The format specifiers `%lu` and `%ld` of the C library are wrong for
`ulong` and `slong` on 64-bit Windows. FLINT's I/O functions
`flint_printf`, `flint_fprintf`, `flint_sprintf`, `flint_snprintf`,
`flint_scanf`, `flint_fscanf` and `flint_sscanf` accept the length
modifier `w` instead: `%wu` and `%wd` print a `ulong` and an `slong` on
every platform (and `%wx` in hexadecimal). With the C library functions,
the macro `WORD_FMT` gives the right length modifier, as in
`printf(WORD_FMT "d", x)`. The `flint_printf` functions also print FLINT
objects directly; see [](flint/printf).

## Platform-dependent behaviour {#platform}

The results of FLINT functions do not depend on the platform, with a few
documented exceptions:

- Functions that take or return words naturally have a word-size range:
  `fmpz_get_si` returns a 64-bit or 32-bit value, `nmod` moduli are
  limited to `FLINT_BITS` bits, and the random functions of
  `ulong_extras` produce word-size values, so random test inputs differ
  between 32-bit and 64-bit platforms.
- Tuning cutoffs, and therefore the choice of algorithms, depend on the
  processor. Results are the same, but timings, and the exact error
  bounds of some ball arithmetic functions, may differ.
- The [](fft_small) module is only available on x86-64 processors with
  AVX2 and on ARM64; elsewhere other algorithms are used.
- `flint_sort` uses `qsort_r` with the GNU C library and a merge sort
  elsewhere, so the order of equal elements may differ.

## Building on Windows {#windows}

On Windows, FLINT can be built with MSYS2 or Cygwin using `configure`,
or natively with CMake and MSVC; see
[](guide/getting-started/building#cmake). With MinGW-w64, FLINT is built
with 64-bit `long long` words (`FLINT_LONG_LONG` is defined in
`flint-config.h`). Programs compiled with MSVC can use the FLINT headers
as long as they are compiled with the same configuration as the library.
