#
#   Copyright (C) 2023 Albin Ahlbäck
#
#   This file is part of FLINT.
#
#   FLINT is free software: you can redistribute it and/or modify it under
#   the terms of the GNU Lesser General Public License (LGPL) as published
#   by the Free Software Foundation; either version 3 of the License, or
#   (at your option) any later version.  See <https://www.gnu.org/licenses/>.
#

include(`config.m4')dnl
dnl
.text

.global	FUNC(flint_mpn_sqr_2)
.p2align	4, 0x90
TYPE(flint_mpn_sqr_2)

FUNC(flint_mpn_sqr_2):
	.cfi_startproc
	mov	0*8(%rsi), %rdx
	xor	%r9d, %r9d
	mulx	1*8(%rsi), %rcx, %r8
	add	%rcx, %rcx
	adc	%r8, %r8
	adc	%r9, %r9
	mulx	%rdx, %r10, %rax
	mov	%r10, 0*8(%rdi)
	add	%rax, %rcx
	mov	1*8(%rsi), %rdx
	mov	%rcx, 1*8(%rdi)
	mulx	%rdx, %r10, %rax
	adc	%r10, %r8
	adc	%r9, %rax
	mov	%r8, 2*8(%rdi)
	mov	%rax, 3*8(%rdi)

	ret
.flint_mpn_sqr_2_end:
SIZE(flint_mpn_sqr_2, .flint_mpn_sqr_2_end)
.cfi_endproc
