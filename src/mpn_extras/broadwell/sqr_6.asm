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

.global	FUNC(flint_mpn_sqr_6)
.p2align	4, 0x90
TYPE(flint_mpn_sqr_6)

FUNC(flint_mpn_sqr_6):
	.cfi_startproc
	mov	0*8(%rsi), %rdx
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13
	push	%r14
	push	%r15
	xor	%r15d, %r15d
	mulx	1*8(%rsi), %rcx, %rbp
	mulx	2*8(%rsi), %r8, %r12
	mulx	3*8(%rsi), %r9, %r13
	mulx	4*8(%rsi), %r10, %r14
	mulx	5*8(%rsi), %r11, %rbx
	adcx	%rbp, %r8
	adcx	%r12, %r9
	adcx	%r13, %r10
	adcx	%r14, %r11
	adcx	%r15, %rbx
	mov	1*8(%rsi), %rdx
	mulx	2*8(%rsi), %r14, %rbp
	adcx	%r14, %r9
	adcx	%rbp, %r10
	mulx	3*8(%rsi), %r14, %rbp
	adox	%r14, %r10
	adox	%rbp, %r11
	mulx	4*8(%rsi), %r14, %rbp
	adcx	%r14, %r11
	adcx	%rbp, %rbx
	mulx	5*8(%rsi), %r14, %rbp
	adox	%r14, %rbx
	adox	%r15, %rbp
	adcx	%r15, %rbp
	mov	2*8(%rsi), %rdx
	mulx	3*8(%rsi), %r14, %r12
	adcx	%r14, %r11
	adcx	%r12, %rbx
	mulx	4*8(%rsi), %r14, %r12
	adox	%r14, %rbx
	adox	%r12, %rbp
	mulx	5*8(%rsi), %r14, %r12
	adcx	%r14, %rbp
	adcx	%r15, %r12
	mov	3*8(%rsi), %rdx
	mulx	4*8(%rsi), %r14, %r13
	adcx	%r14, %rbp
	adcx	%r13, %r12
	mulx	5*8(%rsi), %r14, %r13
	adox	%r14, %r12
	adox	%r15, %r13
	mov	4*8(%rsi), %rdx
	mulx	5*8(%rsi), %rdx, %r14
	adcx	%rdx, %r13
	adc	%r15, %r14
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
	adc	%r14, %r14
	setc	%r15b
	mulx	%rdx, %rdx, %rax
	mov	%rdx, 0*8(%rdi)
	mov	1*8(%rsi), %rdx
	add	%rax, %rcx
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
	adc	%rax, %r13
	mov	5*8(%rsi), %rdx
	mov	%r12, 8*8(%rdi)
	mov	%r13, 9*8(%rdi)
	mulx	%rdx, %rcx, %rax
	adc	%rcx, %r14
	adc	%r15, %rax
	mov	%r14, 10*8(%rdi)
	mov	%rax, 11*8(%rdi)
	pop	%r15
	pop	%r14
	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
.flint_mpn_sqr_6_end:
SIZE(flint_mpn_sqr_6, .flint_mpn_sqr_6_end)
.cfi_endproc
