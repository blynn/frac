#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "cf.h"
#include "test.h"

static void *sqrt2(cf_t cf) {
  cf_put_int(cf, 1);
  while(cf_wait(cf)) {
    cf_put_int(cf, 2);
  }
  return NULL;
}

int main() {
  cf_t x, conv;
  x = cf_new_const(sqrt2);
  conv = cf_new_cf_convergent(x);

  mpz_t p, q;
  mpz_init(p);
  mpz_init(q);

  cf_get(p, conv);
  EXPECT(!mpz_cmp_ui(p, 1));
  cf_get(q, conv);
  EXPECT(!mpz_cmp_ui(q, 1));

  cf_get(p, conv);
  EXPECT(!mpz_cmp_ui(p, 3));
  cf_get(q, conv);
  EXPECT(!mpz_cmp_ui(q, 2));

  cf_get(p, conv);
  EXPECT(!mpz_cmp_ui(p, 7));
  cf_get(q, conv);
  EXPECT(!mpz_cmp_ui(q, 5));

  cf_get(p, conv);
  EXPECT(!mpz_cmp_ui(p, 17));
  cf_get(q, conv);
  EXPECT(!mpz_cmp_ui(q, 12));

  mpz_clear(p);
  mpz_clear(q);
  cf_free(conv);
  cf_free(x);

  x = cf_new_const(sqrt2);
  cf_t mob;
  mpz_t z[4];
  for (int i = 0; i < 4; i++) mpz_init(z[i]);
  // Identity Mobius transformation.
  mpz_set_si(z[0], 1);
  mpz_set_si(z[3], 1);
  mob = cf_new_mobius_to_cf(x, z);
  CF_EXPECT_DEC(mob, "1.4142135623730");
  cf_free(mob);
  cf_free(x);
  x = cf_new_const(sqrt2);
  // (x - 2)/(-3x + 4)
  // Works out to be 1 + sqrt(2).
  mpz_set_si(z[0], 1);
  mpz_set_si(z[1], -2);
  mpz_set_si(z[2], -3);
  mpz_set_si(z[3], 4);
  mob = cf_new_mobius_to_cf(x, z);
  CF_EXPECT_DEC(mob, "2.4142135623730");
  cf_free(x);
  cf_free(mob);
  for (int i = 0; i < 4; i++) mpz_clear(z[i]);

  return 0;
}
