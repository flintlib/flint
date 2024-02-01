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

.global	FUNC(flint_mpn_sqr_3)
.p2align	4, 0x90
TYPE(flint_mpn_sqr_3)

FUNC(flint_mpn_sqr_3):
	.cfi_startproc
	mov	0*8(%rsi), %rdx
	xor	%r11d, %r11d
	mulx	1*8(%rsi), %rcx, %rax
	mulx	2*8(%rsi), %r8, %r9
	mov	1*8(%rsi), %rdx
	add	%rax, %r8
	mulx	2*8(%rsi), %rax, %r10
	adc	%rax, %r9
	adc	%r11, %r10
	mov	0*8(%rsi), %rdx
	add	%rcx, %rcx
	adc	%r8, %r8
	adc	%r9, %r9
	adc	%r10, %r10
	adc	%r11, %r11
	mulx	%rdx, %rdx, %rax
	mov	%rdx, 0*8(%rdi)
	add	%rax, %rcx
	mov	1*8(%rsi), %rdx
	mov	%rcx, 1*8(%rdi)
	mulx	%rdx, %rdx, %rax
	adc	%rdx, %r8
	adc	%rax, %r9
	mov	2*8(%rsi), %rdx
	mov	%r8, 2*8(%rdi)
	mov	%r9, 3*8(%rdi)
	mulx	%rdx, %rdx, %rax
	adc	%rdx, %r10
	adc	%r11, %rax
	mov	%r10, 4*8(%rdi)
	mov	%rax, 5*8(%rdi)

	ret
.flint_mpn_sqr_3_end:
SIZE(flint_mpn_sqr_3, .flint_mpn_sqr_3_end)
.cfi_endproc
