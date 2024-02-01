#
#   Copyright (C) 2024 Albin Ahlbäck
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

.global	FUNC(flint_mpn_mulhigh_6)
.p2align	4, 0x90
TYPE(flint_mpn_mulhigh_6)

FUNC(flint_mpn_mulhigh_6):
	.cfi_startproc
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13

	mov	%rdx, %rcx
	mov	0*8(%rdx), %rdx
	xor	%r9d, %r9d

	mulx	4*8(%rsi), %r8, %rax
	mulx	5*8(%rsi), %r8, %r10
	adcx	%r8, %rax
	adcx	%r9, %r10

	mov	1*8(%rcx), %rdx
	mulx	3*8(%rsi), %r8, %rbx
	mulx	4*8(%rsi), %r8, %r11
	adcx	%rbx, %rax
	adox	%r8, %rax
	adcx	%r11, %r10
	mulx	5*8(%rsi), %r8, %r11
	adox	%r8, %r10
	adcx	%r9, %r11
	adox	%r9, %r11

	mov	2*8(%rcx), %rdx
	mulx	2*8(%rsi), %r8, %rbp
	mulx	3*8(%rsi), %r8, %rbx
	adcx	%rbp, %rax
	adox	%r8, %rax
	adcx	%rbx, %r10
	mulx	4*8(%rsi), %r8, %rbx
	adox	%r8, %r10
	adcx	%rbx, %r11
	mulx	5*8(%rsi), %r8, %rbx
	adox	%r8, %r11
	adcx	%r9, %rbx
	adox	%r9, %rbx

	mov	3*8(%rcx), %rdx
	mulx	1*8(%rsi), %r8, %r12
	mulx	2*8(%rsi), %r8, %rbp
	adcx	%r12, %rax
	adox	%r8, %rax
	adcx	%rbp, %r10
	mulx	3*8(%rsi), %r8, %rbp
	adox	%r8, %r10
	adcx	%rbp, %r11
	mulx	4*8(%rsi), %r8, %rbp
	adox	%r8, %r11
	adcx	%rbp, %rbx
	mulx	5*8(%rsi), %r8, %rbp
	adox	%r8, %rbx
	adcx	%r9, %rbp
	adox	%r9, %rbp

	mov	4*8(%rcx), %rdx
	mulx	0*8(%rsi), %r8, %r13
	mulx	1*8(%rsi), %r8, %r12
	adcx	%r13, %rax
	adox	%r8, %rax
	adcx	%r12, %r10
	mulx	2*8(%rsi), %r8, %r12
	adox	%r8, %r10
	adcx	%r12, %r11
	mulx	3*8(%rsi), %r8, %r12
	adox	%r8, %r11
	adcx	%r12, %rbx
	mulx	4*8(%rsi), %r8, %r12
	adox	%r8, %rbx
	adcx	%r12, %rbp
	mulx	5*8(%rsi), %r8, %r12
	adox	%r8, %rbp
	adcx	%r9, %r12
	adox	%r9, %r12

	mov	5*8(%rcx), %rdx
	mulx	0*8(%rsi), %r8, %r13
	adcx	%r8, %rax
	adcx	%r13, %r10
	mulx	1*8(%rsi), %r8, %r13
	adox	%r8, %r10
	adox	%r13, %r11
	mulx	2*8(%rsi), %r8, %r13
	adcx	%r8, %r11
	adcx	%r13, %rbx
	mulx	3*8(%rsi), %r8, %r13
	adox	%r8, %rbx
	adox	%r13, %rbp
	mulx	4*8(%rsi), %r8, %r13
	adcx	%r8, %rbp
	adcx	%r13, %r12
	mulx	5*8(%rsi), %r8, %r13
	adox	%r8, %r12
	adcx	%r9, %r13
	adox	%r9, %r13

	mov	%r10, 0*8(%rdi)
	mov	%r11, 1*8(%rdi)
	mov	%rbx, 2*8(%rdi)
	mov	%rbp, 3*8(%rdi)
	mov	%r12, 4*8(%rdi)
	mov	%r13, 5*8(%rdi)

	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
.flint_mpn_mulhigh_6_end:
SIZE(flint_mpn_mulhigh_6, .flint_mpn_mulhigh_6_end)
.cfi_endproc
