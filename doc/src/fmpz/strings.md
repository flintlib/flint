---
title: Strings, printing and binary input and output
functions: fmpz_get_str fmpz_set_str fmpz_fprint fmpz_fread fmpz_out_raw
---

`fmpz_get_str` and `fmpz_set_str` convert between integers and strings in
any base from 2 to 62, `fmpz_print` and `fmpz_fprint` print in decimal,
and `fmpz_read` and `fmpz_fread` parse decimal input from a stream.
`fmpz_out_raw` and `fmpz_inp_raw` store integers in GMP's portable binary
format.

## Formats {#formats}

Digits are `0`–`9` for bases up to 10, then lower-case letters up to base
36 (`fmpz_set_str` accepts both cases there), and for bases 37 to 62
`0`–`9`, `A`–`Z`, `a`–`z`, as in GMP. A negative number starts with a
minus sign; there is no plus sign and no prefix such as `0x`.

`fmpz_get_str(str, b, f)` writes into `str` if it is not `NULL`, which
must have room for `fmpz_sizeinbase(f, b) + 2` characters (the digits,
a sign and the terminating zero); `fmpz_sizeinbase` is exact for every
base. With `str = NULL` the string is allocated with `flint_malloc` and
must be freed with `flint_free`. The function returns the string.

`fmpz_set_str` accepts white space before the number, between the minus
sign and the digits, and after the number; white space between digits is
ignored as well. It returns 0 on success and $-1$ if the string is not a
valid integer, in which case the value of the output is undefined.

`fmpz_fprint` writes the decimal representation without a newline and
returns the number of characters written. `fmpz_fread` skips leading
white space and reads an optional minus sign and decimal digits; it
returns a positive value on success. The raw format of `fmpz_out_raw` is
a four-byte big-endian signed byte count followed by the magnitude in
big-endian byte order; it is portable between machines and compatible
with `mpz_out_raw` and `mpz_inp_raw`.

## Algorithms {#algorithms}

In base 10, inline values are converted with word arithmetic (other
bases use GMP's `mpz_get_str` and `mpz_set_str` for all sizes). For
large values:

- `fmpz_get_str` in base 10 uses FLINT's radix conversion
  (`flint_mpn_get_str`, built on the [](radix) module) for up to 40 and
  from 1400 limbs on in builds with `fft_small`, and GMP's
  `mpz_get_str` in between, for other bases and in builds without
  `fft_small`;
- `fmpz_set_str` in base 10 uses GMP's `mpn_set_str` below 24000 digits
  and above that a divide-and-conquer conversion, which splits the string
  in halves, converts them (in parallel if threads are available) and
  combines them as $h \cdot 10^k + l$ with precomputed powers of 5 (since
  $10^k = 5^k 2^k$); other bases, and strings with characters other
  than digits after the sign (such as interior white space), go to
  `mpz_set_str`.

Both large-number algorithms run in $O(M(n) \log n)$ for $n$ digits,
where $M(n)$ is the cost of multiplication;
[](algorithms/integers/radix-conversion) describes them.

## Example {#example}

```c example name=fmpz_strings
#include <stdio.h>
#include <flint/flint.h>
#include <flint/fmpz.h>

int main(void)
{
    fmpz_t x, y;
    char * s;
    char buf[100];
    int bases[4] = { 2, 16, 36, 62 }, i;
    FILE * f;

    fmpz_init(x);
    fmpz_init(y);

    fmpz_set_str(x, "-123456789012345678901234567890", 10);
    for (i = 0; i < 4; i++)
    {
        s = fmpz_get_str(NULL, bases[i], x);    /* allocated by FLINT */
        flint_printf("base %2d: %s\n", bases[i], s);
        fmpz_set_str(y, s, bases[i]);          /* and back */
        flint_free(s);
        if (!fmpz_equal(x, y))
            flint_printf("round trip failed\n");
    }

    /* a caller-supplied buffer needs sizeinbase + 2 characters */
    fmpz_one_2exp(x, 100);
    flint_printf("2^100 needs %zu digits: %s\n",
                 fmpz_sizeinbase(x, 10), fmpz_get_str(buf, 10, x));

    /* white space and invalid input */
    flint_printf("\"  -42 \" -> %d, ", fmpz_set_str(y, "  -42 ", 10));
    fmpz_print(y); flint_printf("\n");
    flint_printf("\"12a\" -> %d\n", fmpz_set_str(y, "12a", 10));

    /* binary round trip through a file in GMP's raw format */
    f = tmpfile();
    if (f != NULL)
    {
        flint_printf("raw: %zu bytes written, ", fmpz_out_raw(f, x));
        rewind(f);
        flint_printf("%zu read, ", fmpz_inp_raw(y, f));
        flint_printf("equal: %d\n", fmpz_equal(x, y));
        fclose(f);
    }

    fmpz_clear(x);
    fmpz_clear(y);
    return 0;
}
```
```output
base  2: -1100011101110100100001111111101101100001101110011111000001110111001001110001111110000101011010010
base 16: -18ee90ff6c373e0ee4e3f0ad2
base 36: -byw97um9s91dlz68tsi
base 62: -2aYls9bkamJJSwhr0
2^100 needs 31 digits: 1267650600228229401496703205376
"  -42 " -> 0, -42
"12a" -> -1
raw: 17 bytes written, 17 read, equal: 1
```

The raw representation of $2^{100}$ consists of the four-byte count and
13 bytes of magnitude.
