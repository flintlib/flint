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

.global	FUNC(flint_mpn_sqr_5)
.p2align	4, 0x90
TYPE(flint_mpn_sqr_5)

FUNC(flint_mpn_sqr_5):
	.cfi_startproc
	mov	0*8(%rsi), %rdx
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13
	xor	%r13d, %r13d
	mulx	1*8(%rsi), %rcx, %rbx
	mulx	2*8(%rsi), %r8, %rbp
	mulx	3*8(%rsi), %r9, %r12
	mulx	4*8(%rsi), %r10, %r11
	adcx	%rbx, %r8
	adcx	%rbp, %r9
	adcx	%r12, %r10
	adcx	%r13, %r11
	mov	1*8(%rsi), %rdx
	mulx	2*8(%rsi), %r12, %rbx
	adcx	%r12, %r9
	adcx	%rbx, %r10
	mulx	3*8(%rsi), %r12, %rbx
	adox	%r12, %r10
	adox	%rbx, %r11
	mulx	4*8(%rsi), %r12, %rbx
	adcx	%r12, %r11
	adcx	%r13, %rbx
	mov	2*8(%rsi), %rdx
	mulx	3*8(%rsi), %r12, %rbp
	adcx	%r12, %r11
	adcx	%rbp, %rbx
	mulx	4*8(%rsi), %r12, %rbp
	adox	%r12, %rbx
	adox	%r13, %rbp
	mov	3*8(%rsi), %rdx
	mulx	4*8(%rsi), %rdx, %r12
	adcx	%rdx, %rbp
	adc	%r13, %r12
	mov	0*8(%rsi), %rdx
	add	%rcx, %rcx
	adc	%r8, %r8
	adc	%r9, %r9
	adc	%r10, %r10
	adc	%r11, %r11
	adc	%rbx, %rbx
	adc	%rbp, %rbp
	adc	%r12, %r12
	adc	%r13, %r13
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
	mulx	%rdx, %rcx, %rax
	adc	%rcx, %r10
	adc	%rax, %r11
	mov	3*8(%rsi), %rdx
	mov	%r10, 4*8(%rdi)
	mov	%r11, 5*8(%rdi)
	mulx	%rdx, %rcx, %rax
	adc	%rcx, %rbx
	adc	%rax, %rbp
	mov	4*8(%rsi), %rdx
	mov	%rbx, 6*8(%rdi)
	mov	%rbp, 7*8(%rdi)
	mulx	%rdx, %rcx, %rax
	adc	%rcx, %r12
	adc	%r13, %rax
	mov	%r12, 8*8(%rdi)
	mov	%rax, 9*8(%rdi)
	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
.flint_mpn_sqr_5_end:
SIZE(flint_mpn_sqr_5, .flint_mpn_sqr_5_end)
.cfi_endproc
