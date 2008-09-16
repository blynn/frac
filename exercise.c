// Compute 37 digits of sqrt(2).

#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {
  // 4 byte precision gives 18 digits,
  // 8 byte precision gives 37 digits.
  int n = 37;
  if (argc > 1) {
    n = atoi(argv[1]);
    if (n <= 0) n = 37;
  }
  unsigned long long p0, p1, q0, q1, r0, r1;

  // sqrt(2) = [1; 2, 2, ...]
  // Handle first term here.
  p0 = p1 = q1 = 1;
  q0 = 0;

  while (n > 0) {
    // Continued fraction recurrence relations.
    r0 = p0 + (p1 << 1);
    p0 = p1;
    p1 = r0;

    r0 = q0 + (q1 << 1);
    q0 = q1;
    q1 = r0;

    // digit <- floor(p0 / q0)
    // r0 = q0 * (digit + 1)
    // r1 = q1 * (digit + 1)
    int digit = 0;
    r0 = q0;
    r1 = q1;
    while (digit < 9 && r0 <= p0) {
      r0 += q0;
      r1 += q1;
      digit++;
    }

    if (r1 > p1) {
      r1 -= q1;
      if (r1 <= p1) {
	n--;
	putchar(digit + '0');
	p1 -= r1;
	r0 -= q0;
	p0 -= r0;

	// p0 <- 10 * p0, p1 <- 10 * p1
	p0 = ((p0 << 2) + p0) << 1;
	p1 = ((p1 << 2) + p1) << 1;
      }
    }
  }
  putchar('\n');
  return 0;
}
