// Test quadratic algorithm.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "cf.h"
#include "test.h"

int main() {
  // I used bc to get the reference values.
  // Check e + pi
  cf_t e = cf_new_e();
  cf_t pi = cf_new_pi();

  cf_t b = cf_new_add(e, pi);

  CF_EXPECT_DEC(b, "5.8598744820488384738");

  cf_free(b);
  cf_free(e);
  cf_free(pi);

  // Check e * pi
  e = cf_new_e();
  pi = cf_new_pi();

  b = cf_new_mul(e, pi);

  CF_EXPECT_DEC(b, "8.53973422267356706546");

  cf_free(b);
  cf_free(e);
  cf_free(pi);

  // Check 2 sin 1 cos 1 = sin 2.
  cf_t s1, c1;
  s1 = cf_new_sin1();
  c1 = cf_new_cos1();
  mpz_t a[8];
  mpz8_init(a);
  mpz8_set_int(a,
      2, 0, 0, 0,
      0, 0, 0, 1);
  b = cf_new_bihom(s1, c1, a);

  CF_EXPECT_DEC(b, "0.9092974268256816953");

  cf_free(b);
  cf_free(c1);
  cf_free(s1);

  // Check 2 (cos 1)^2 - 1 = cos 2.
  c1 = cf_new_cos1();
  cf_t t[2];
  cf_tee(t, c1);
  mpz8_set_int(a,
      2, 0, 0, -1,
      0, 0, 0, 1);
  b = cf_new_bihom(t[0], t[1], a);

  CF_EXPECT_DEC(b, "-0.41614683654714238699");

  cf_free(t[0]);
  cf_free(t[1]);
  cf_free(b);
  cf_free(c1);

  mpz8_clear(a);
  return 0;
}
