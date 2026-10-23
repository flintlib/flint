---
title: Proving primality with certificates
functions: ecpp_prove ecpp_verify ecpp_is_prime
---

`ecpp_prove` proves that a probable prime $n$ is prime and returns a
certificate; `ecpp_verify` checks a certificate. The verification only
needs a few scalar multiplications on elliptic curves modulo the numbers
of the chain and is typically one to two orders of magnitude faster than
the proof. Certificates can be exchanged with PARI/GP (`primecert` and
`primecertisvalid`) through the text formats.

::: summary algorithms/number-theory/primality-proving :::

## When to use ECPP {#when}

For a yes/no answer, call `fmpz_is_prime`, which decides small inputs and
inputs with easily factored $n \pm 1$ directly and otherwise runs APRCL
or ECPP, whichever is expected to be faster: ECPP from 250 bits with one
thread, 550 with two, 1000 with three, 2000 with four and 8000 with up to
eight threads. APRCL parallelises better (its Jacobi sum tests are
independent) while the ECPP chain is largely serial, so the crossover
moves up with the number of threads. Use `ecpp_prove` directly when a
certificate is wanted.

For example, proving $10^{100} + 267$ took 17 ms with ECPP and 29 ms with
APRCL on one core of the machine used for these timings, and
$2^{521} - 1$ (157 digits) 84 ms and 149 ms.

## How the proof proceeds {#steps}

Each step replaces $n$ by a probable prime $q$ roughly $\log_2 |D|$ bits
smaller: it finds a discriminant $D$ for which $4n = t^2 + |D| v^2$ has a
solution and one of the possible curve orders $m = n + 1 \pm t$ is a small
number times a probable prime $q > (n^{1/4} + 1)^2$, constructs a curve
with $m$ points by complex multiplication (a root of the Hilbert class
polynomial $H_D$ modulo $n$, found through the class field tower or the
genus field), and finds a point $P$ with $mP = O$ and $(m/q) P \ne O$.
The search for the next step is organised as in fastECPP: discriminants
smooth over a set of small primes, batch trial division of all candidate
orders by a remainder tree, and the candidate with the best estimated cost
per bit gained. The chain ends when $q < 2^{64}$.

With threads, the square roots, the Cornacchia steps, the probable prime
tests and the reduction of the primorial are distributed, the two twists
of a curve are tried at once, and the realisation of a step (class
polynomial, root, curve and point) runs on a pool thread while the search
for the next step proceeds.

`ecpp_set_verbose(1)` reports the steps as they are found. The example
program `factor_integer` (in `examples/`) proves the prime factors it
finds with ECPP when given `-certify`, and prints the certificates in the
PARI/GP format with `-pari`.

## Prior art {#prior-art}

The organisation of the search follows fastECPP of Franke, Kleinjung,
Morain and Wirth [@FKMW2004] and Morain [@Morain2007], in the form of the
implementations by Jared Asuncion in PARI/GP (`primecert`, [@PARI]) and by
Andreas Enge in the CM library [@EngeCM; @Enge2024]: discriminants smooth
over a set of small primes, square roots of those primes computed once
per $n$, batch factoring of the curve orders with a remainder tree, and
the choice of the trial division bound and of the minimum gain per step.
The rounds of the search, which extend the set of primes until a few
prime cofactors are expected, the trial division bound fixed for the
whole chain, the cache-blocked class number sieve and the class field
tower with the Hecke representation of its coefficients are as in CM,
which also provided the reference for the normalisations of Weber's
function as a class invariant [@Schertz2002; @YuiZagier1997]. The class
field tower is the decomposition of Enge and Morain [@EngeMorain2003];
the factor of $H_D$ over the genus field and the descent by radicals and
by Kummer theory through the tower are original to this implementation.

## Example {#example}

Prove that $10^{100} + 267$ is prime, check the certificate, and read it
back from the PARI/GP format:

```c example name=ecpp_prove_verify
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/ecpp.h>

int main(void)
{
    fmpz_t n;
    ecpp_cert_t cert, cert2;
    char * s;
    int ok;

    fmpz_init(n);
    ecpp_cert_init(cert);
    ecpp_cert_init(cert2);

    fmpz_ui_pow_ui(n, 10, 100);
    fmpz_add_ui(n, n, 267);

    flint_printf("proved: %d\n", ecpp_prove(cert, n));
    flint_printf("verified: %d\n", ecpp_verify(cert, n));
    flint_printf("first step: n = %{fmpz}\n", cert->steps[0].n);
    flint_printf("last q below 2^64: %d\n",
        fmpz_bits(cert->steps[cert->num - 1].q) <= 64);

    s = ecpp_cert_get_str(cert, ECPP_CERT_FORMAT_PARI);
    ok = ecpp_cert_set_str(cert2, s);
    flint_printf("parsed: %d, verified: %d\n", ok, ecpp_verify(cert2, n));

    /* the certificate does not prove anything about another number */
    fmpz_add_ui(n, n, 2);
    flint_printf("other n: %d\n", ecpp_verify(cert, n));

    flint_free(s);
    ecpp_cert_clear(cert);
    ecpp_cert_clear(cert2);
    fmpz_clear(n);
    return 0;
}
```
```output
proved: 1
verified: 1
first step: n = 10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000267
last q below 2^64: 1
parsed: 1, verified: 1
other n: 0
```
