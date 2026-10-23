---
title: Fast Fourier transforms
short: FFT
summary: The FFTs behind fast multiplication of integers and polynomials: Schönhage–Strassen, small-prime FFTs in floating point, and generic DFTs over any ring.
order: 8
---

FLINT has three FFT implementations, each suited to different coefficient
rings and sizes. The sections below describe the transforms, their truncated
and cache-friendly variants, and how they are combined with the Chinese
remainder theorem to multiply large integers and polynomials.
