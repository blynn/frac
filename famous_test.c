// Test famous continued fraction expansions are correct.
// as well as expansions based on Taylor series.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "cf.h"
#include "test.h"

int main() {
  CF_NEW_EXPECT_DEC(cf_new_sqrt2, "1.41421356237309504880");
  CF_NEW_EXPECT_DEC(cf_new_e, "2.71828182845904523536");
  CF_NEW_EXPECT_DEC(cf_new_pi, "3.1415926535897932384");
  CF_NEW_EXPECT_DEC(cf_new_tan1, "1.5574077246549022305");
  CF_NEW_EXPECT_DEC(cf_new_sin1, "0.8414709848078965066");
  CF_NEW_EXPECT_DEC(cf_new_cos1, "0.5403023058681397174");

  mpz_t z;
  mpz_init(z);
  cf_t x;

  mpz_set_si(z, 3);
  x = cf_new_tanh(z);
  CF_EXPECT_DEC(x, "0.99505475368673045133188018525548847509781385470028");
  cf_free(x);

  mpz_set_ui(z, 2);
  x = cf_new_epow(z);
  CF_EXPECT_DEC(x, "7.38905609893065022723042746057500781318031557055184");
  cf_free(x);

  mpz_clear(z);
  return 0;
}
