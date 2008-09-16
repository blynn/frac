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
  mpz_t addarray[8];
  for (int i = 0; i < 8; i++) {
    mpz_init(addarray[i]);
  }
  mpz_set_ui(addarray[1], 1);
  mpz_set_ui(addarray[2], 1);
  mpz_set_ui(addarray[7], 1);

  cf_t b = cf_new_bihom(e, pi, addarray);

  CF_EXPECT_DEC(b, "5.8598744820488384738");

  cf_free(b);
  cf_free(e);
  cf_free(pi);

  // Check 2 sin 1 cos 2 = sin 2 = 0.909...
  cf_t s1, c1;
  s1 = cf_new_sin1();
  c1 = cf_new_cos1();
  mpz_set_ui(addarray[0], 2);
  mpz_set_ui(addarray[1], 0);
  mpz_set_ui(addarray[2], 0);
  mpz_set_si(addarray[7], 1);
  b = cf_new_bihom(s1, c1, addarray);

  CF_EXPECT_DEC(b, "0.9092974268256816953");

  cf_free(b);
  cf_free(c1);
  cf_free(s1);

  // Check 2 (cos 1)^2 - 1 = cos 2 =
  s1 = cf_new_cos1();  // TODO: Implement tee, use that instead.
  c1 = cf_new_cos1();
  mpz_set_si(addarray[0], 2);
  mpz_set_si(addarray[3], -1);
  mpz_set_si(addarray[7], 1);
  b = cf_new_bihom(s1, c1, addarray);

  CF_EXPECT_DEC(b, "-0.41614683654714238699");

  cf_free(b);
  cf_free(c1);
  cf_free(s1);

  for (int i = 0; i < 8; i++) {
    mpz_clear(addarray[i]);
  }
  return 0;
}
