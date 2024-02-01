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

.global	FUNC(flint_mpn_sqr_7)
.p2align	4, 0x90
TYPE(flint_mpn_sqr_7)

FUNC(flint_mpn_sqr_7):
	.cfi_startproc
	mov	0*8(%rsi), %rdx
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13
	push	%r14
	push	%r15
	xor	%r10d, %r10d
	mulx	1*8(%rsi), %rcx, %r11
	mulx	2*8(%rsi), %r8, %rbx
	mulx	3*8(%rsi), %r9, %r14
	adcx	%r11, %r8
	adcx	%rbx, %r9
	mulx	4*8(%rsi), %rax, %r12
	mulx	5*8(%rsi), %r11, %r13
	mulx	6*8(%rsi), %rbx, %rbp
	adcx	%r14, %rax
	adcx	%r12, %r11
	mov	1*8(%rsi), %rdx
	adcx	%r13, %rbx
	adcx	%r10, %rbp
	mulx	2*8(%rsi), %r14, %r12
	mulx	3*8(%rsi), %r15, %r13
	adcx	%r14, %r9
	adcx	%r12, %rax
	adox	%r15, %rax
	adox	%r13, %r11
	mulx	4*8(%rsi), %r14, %r12
	mulx	5*8(%rsi), %r15, %r13
	adcx	%r14, %r11
	adcx	%r12, %rbx
	mulx	6*8(%rsi), %r14, %r12
	adox	%r15, %rbx
	adox	%r13, %rbp
	mov	2*8(%rsi), %rdx
	adcx	%r14, %rbp
	adox	%r10, %r12
	adcx	%r10, %r12
	mulx	3*8(%rsi), %r15, %r13
	adox	%r15, %r11
	adox	%r13, %rbx
	mulx	4*8(%rsi), %r15, %r13
	adcx	%r15, %rbx
	adcx	%r13, %rbp
	mulx	5*8(%rsi), %r15, %r13
	adox	%r15, %rbp
	adox	%r13, %r12
	mulx	6*8(%rsi), %r15, %r13
	adcx	%r15, %r12
	adox	%r10, %r13
	mov	3*8(%rsi), %rdx
	adcx	%r10, %r13
	mulx	4*8(%rsi), %r15, %r14
	adcx	%r15, %rbp
	adcx	%r14, %r12
	mulx	5*8(%rsi), %r15, %r14
	adox	%r15, %r12
	adox	%r14, %r13
	mulx	6*8(%rsi), %r15, %r14
	mov	4*8(%rsi), %rdx
	adcx	%r15, %r13
	adcx	%r10, %r14
	mulx	5*8(%rsi), %r10, %r15
	adcx	%r10, %r13
	adcx	%r15, %r14
	mulx	6*8(%rsi), %r10, %r15
	mov	$0, %edx
	adox	%r10, %r14
	adox	%rdx, %r15
	mov	5*8(%rsi), %rdx
	mulx	6*8(%rsi), %rdx, %r10
	adcx	%rdx, %r15
	adc	$0, %r10
	test	%al, %al
	mov	0*8(%rsi), %rdx
	adcx	%rcx, %rcx
	adcx	%r8, %r8
	adcx	%r9, %r9
	adcx	%rax, %rax
	adcx	%r11, %r11
	adcx	%rbx, %rbx
	push	%rax
	adcx	%rbp, %rbp
	adcx	%r12, %r12
	adcx	%r13, %r13
	adcx	%r14, %r14
	adcx	%r15, %r15
	adcx	%r10, %r10
	mulx	%rdx, %rdx, %rax
	mov	%rdx, 0*8(%rdi)
	adox	%rax, %rcx
	mov	1*8(%rsi), %rdx
	mov	%rcx, 1*8(%rdi)
	mulx	%rdx, %rdx, %rax
	adox	%rdx, %r8
	adox	%rax, %r9
	pop	%rcx
	mov	2*8(%rsi), %rdx
	mov	%r8, 2*8(%rdi)
	mov	%r9, 3*8(%rdi)
	mulx	%rdx, %rdx, %rax
	adox	%rdx, %rcx
	adox	%rax, %r11
	mov	3*8(%rsi), %rdx
	mov	%rcx, 4*8(%rdi)
	mov	%r11, 5*8(%rdi)
	mulx	%rdx, %r8, %rax
	adox	%r8, %rbx
	adox	%rax, %rbp
	mov	4*8(%rsi), %rdx
	mov	%rbx, 6*8(%rdi)
	mov	%rbp, 7*8(%rdi)
	mulx	%rdx, %r8, %rax
	adox	%r8, %r12
	adox	%rax, %r13
	mov	5*8(%rsi), %rdx
	mov	%r12, 8*8(%rdi)
	mov	%r13, 9*8(%rdi)
	mulx	%rdx, %r8, %rax
	adox	%r8, %r14
	adox	%rax, %r15
	mov	$0, %ecx
	mov	6*8(%rsi), %rdx
	mov	%r14, 10*8(%rdi)
	mov	%r15, 11*8(%rdi)
	mulx	%rdx, %r8, %rax
	adox	%r8, %r10
	adox	%rcx, %rax
	adcx	%rcx, %rax
	mov	%r10, 12*8(%rdi)
	mov	%rax, 13*8(%rdi)
	pop	%r15
	pop	%r14
	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
.flint_mpn_sqr_7_end:
SIZE(flint_mpn_sqr_7, .flint_mpn_sqr_7_end)
.cfi_endproc
