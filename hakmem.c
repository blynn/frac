// Compute the HAKMEM constant.
//   (sqrt(3/pi^2 + e)/(tanh(sqrt(5))-sin(69))
//
// Confirm with:
//  $ echo "pi=4*a(1)
//          x=2*sqrt(5)
//          (sqrt(3/pi^2+e(1)))/((e(x)-1)/(e(x)+1)-s(69))" | bc -l
//
// Almost exclusively uses techniques described by Gosper. Differences:
//
//  - Taylor series for cos(1) to find its continued fraction expansion.
//  - Binary search instead of Newton's method when finding integer part of
//  solution of quadratic.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "cf.h"

void *tanhsqrt5denom(cf_t cf) {
  int odd = 1;
  while(cf_wait(cf)) {
    cf_put_int(cf, odd);
    odd += 2;
    cf_put_int(cf, 5);
  }
  return NULL;
}

void cf_dump(cf_t cf, int n) {
  mpz_t z;
  mpz_init(z);
  cf_t dec = cf_new_cf_to_decimal(cf);
  int i;
  for (i = 0; i <= n; i++) {
    cf_get(z, dec);
    gmp_printf("%Zd", z);
    if (!(i % 5)) putchar(' ');
    if (!(i % 50)) putchar('\n');
  }
  if (n % 50) putchar('\n');
  cf_free(dec);
  mpz_clear(z);
}

int main(int argc, char **argv) {
  int n = 100;
  if (argc > 1) {
    n = atoi(argv[1]);
    if (n <= 0) n = 100;
  }
  cf_t c[7];
  cf_t t[6][2];
  mpz_t b[8];
  int i;

  mpz8_init(b);
  mpz8_set_int(b,
      2, 0, 0, -1,
      0, 0, 0, 1);

  c[0] = cf_new_cos1();
  for (i = 1; i < 7; i++) {
    cf_tee(t[i - 1], c[i - 1]);
    c[i] = cf_new_bihom(t[i - 1][0], t[i - 1][1], b);
  }
  // c[6] = cos 64

  cf_t c1 = cf_new_cos1();
  cf_t ct0[2], ct1[2], ct2[2];
  cf_tee(ct1, c1);
  cf_tee(ct2, ct1[1]);
  cf_t c2 = cf_new_mul(ct2[0], ct2[1]);
  cf_tee(ct0, c2);
  mpz8_set_int(b,
      16, -20, 0, 5,  // cos 5n (5th Chebyshev polynomial)
      0, 0, 0, 1);
  cf_t c3 = cf_new_bihom(ct0[0], ct0[1], b);
  cf_t c5 = cf_new_mul(ct1[0], c3);
  // c5 = cos 5

  cf_t c64t[2], c64t1[2];
  cf_tee(c64t, c[6]);
  cf_tee(c64t1, c64t[1]);
  mpz8_set_int(b,
      -1, 0, 0, 1,
       0, 0, 0, 1);
  cf_t s1 = cf_new_bihom(c64t1[0], c64t1[1], b);
  cf_t s64 = cf_new_sqrt(s1);
  // s64 = sin 64, c64t[0] = untouched cos 64

  cf_t c5t[2], c5t1[2];
  cf_tee(c5t, c5);
  cf_tee(c5t1, c5t[1]);
  cf_t s3 = cf_new_bihom(c5t1[0], c5t1[1], b);
  cf_t s5 = cf_new_sqrt(s3);
  // s5 = sin 5, c5t[0] = untouched cos 5

  cf_t cf0 = cf_new_mul(s5, c64t[0]);
  cf_t cf1 = cf_new_mul(s64, c5t[0]);
  cf_t s69 = cf_new_sub(cf0, cf1); // TODO: Respect signs.
				   // This is supposed to be an addition,
				   // and the result should be negative.
  // s69 = sin 69

  cf_t sqrt5 = cf_new_sqrt5();
  cf_t ts5d = cf_new_const_nonregular(tanhsqrt5denom);
  cf_t ts5 = cf_new_div(sqrt5, ts5d);
  
  cf_t den = cf_new_add(ts5, s69); // TODO: This should be a subtract, but
                                   // again, we don't handle signs properly.

  cf_t e = cf_new_e();
  cf_t pi = cf_new_pi();
  cf_t pitee[2];
  cf_tee(pitee, pi);
  mpz8_set_int(b,
      0, 0, 0, 3,
      1, 0, 0, 0);
  // tops = Three-Over-Pi-Squared.
  cf_t tops = cf_new_bihom(pitee[0], pitee[1], b);

  mpz8_set_add(b);
  cf_t sum = cf_new_bihom(tops, e, b);
  cf_t num = cf_new_sqrt(sum);
  cf_t hakmem_constant = cf_new_div(num, den);
  cf_dump(hakmem_constant, n);

  mpz8_clear(b);
  return 0;
}
