---
title: Security
summary: FLINT is a library for mathematical research, education and experimentation; why it must not be used to implement cryptographic protocols, and why untrusted input should only be processed in a sandbox.
order: 5
---

FLINT is written for mathematical research, education and
experimentation, where the user controls the input and wants correct
results as fast as possible. It is not designed to be safe in an
adversarial setting.

## Cryptography {#cryptography}

Cryptographic implementations must resist attackers who observe or
influence the computation. FLINT makes no attempt to:

- *Run in constant time.* Running times and memory access patterns
  depend on the data. For example, `n_powmod2_ui_preinv` branches on
  each bit of the exponent, `fmpz_powm` calls GMP's `mpz_powm` and not
  its side-channel resistant `mpz_powm_sec`, GCD and division algorithms
  take data-dependent paths, and many functions choose an algorithm
  according to the sizes of their operands. Such timing differences can
  reveal secret exponents and keys.
- *Protect secrets in memory.* Memory is released with `flint_free`,
  which does not overwrite it, and the pool of `fmpz` integers keeps
  freed integers, with their digits, for reuse.
- *Produce unpredictable random numbers.* `flint_rand_t` is a pair of
  linear congruential generators with fixed parameters, which are
  adequate for testing but predictable, and `flint_rand_init` always
  starts from the same seed. Random primes from `n_randprime` or `fmpz_randprime`, random
  residues and every other random object derived from it are therefore
  predictable. The generator exists to produce test inputs and to make
  randomised algorithms reproducible.

Use a dedicated cryptographic library for cryptographic protocols.
Experimenting with cryptographic algorithms (factoring, discrete logarithms,
lattice reduction, prototypes of schemes) is a legitimate use
of FLINT; deploying them is not.

## Untrusted input {#untrusted}

Input from an untrusted source, such as a network service that accepts
expressions, a file uploaded by a user, or a computation requested
through a web interface, must not be passed to FLINT in a process that
has access to anything valuable. The reasons are the topics of the
previous pages:

- *No input validation.* Most functions assume valid input, and the
  preconditions that are checked at all are mostly checked only in
  assertion builds; see [](guide/correctness/input-validation). Invalid
  input can corrupt memory.
- *No overflow checks.* Sizes, shifts, exponents and precisions are not
  checked against overflow; see [](guide/correctness/limits). A small
  input can request an impossible amount of memory, or wrap around in a
  size computation: `fmpz_ui_pow_ui(x, 3, e)` with a 19-digit `e` aborts
  the process when the allocation fails.
- *Unbounded running time.* The cost of many operations is not bounded
  by a function of the size of the input that is reasonable to wait
  for: the factorisation of a 500-digit integer, or a zero test of a
  complicated exact expression.
- *Errors terminate the process.* By default, an exception raised with
  `flint_throw`, and a failed allocation, abort the process.
  A server that calls FLINT directly can be stopped by a single bad
  request.
- *Parsers are not hardened.* The string functions trust the sizes in
  their input. `fmpz_poly_set_str` and `nmod_poly_set_str` allocate the
  number of coefficients announced at the start of the string before
  reading any of them, so `"1000000000000  1"` requests 8 TB of memory.
  `fmpz_poly_set_str` reads past the end of a string announcing more
  coefficients than it contains, and `nmod_poly_set_str` crashes on a
  negative length.
  `nmod_poly_set_str` also accepts coefficients that are not reduced
  modulo $n$, which produces an invalid polynomial.

The safe way to compute with untrusted input is to do it in a separate
process with limited privileges and resource limits: a container or
another sandbox, a limit on the address space (so that allocations fail
instead of exhausting the machine) and on the CPU time, and a supervisor
that treats a crash or a timeout as a failed request. Sizes, degrees,
exponents and precisions should be checked against explicit bounds
before they are passed to FLINT. Language bindings and computer algebra
systems that expose FLINT to their users usually follow this model
already: the user controls the process, and can only harm themselves.
