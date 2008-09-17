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

// Converges extremely slowly.
static void *slow_pi(cf_t cf) {
  mpz_t num, denom, t;
  mpz_init(num);
  mpz_init(denom);
  mpz_init(t);

  mpz_set_ui(denom, 3);
  cf_put(cf, denom);

  mpz_set_ui(num, 1);
  cf_put(cf, num);
  mpz_set_ui(t, 8);

  mpz_set_ui(denom, 6);
  while(cf_wait(cf)) {
    cf_put(cf, denom);
    mpz_add(num, num, t);
    cf_put(cf, num);
    mpz_add_ui(t, t, 8);
  }

  mpz_clear(num);
  mpz_clear(denom);
  mpz_clear(t);
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

  mpz_t z[4];
  for (int i = 0; i < 4; i++) mpz_init(z[i]);
  // Identity Mobius transformation.
  mpz_set_si(z[0], 1);
  mpz_set_si(z[3], 1);
  mpz_t digit;
  mpz_init(digit);

  x = cf_new_const(slow_pi);
  conv = cf_new_nonregular_mobius_to_decimal(x, z);
  cf_get(digit, conv);
  EXPECT(!mpz_cmp_ui(digit, 3));
  cf_get(digit, conv);
  EXPECT(!mpz_cmp_ui(digit, 1));
  cf_get(digit, conv);
  EXPECT(!mpz_cmp_ui(digit, 4));
  cf_get(digit, conv);
  EXPECT(!mpz_cmp_ui(digit, 1));
  cf_get(digit, conv);
  EXPECT(!mpz_cmp_ui(digit, 5));
  cf_get(digit, conv);
  EXPECT(!mpz_cmp_ui(digit, 9));
  cf_get(digit, conv);
  EXPECT(!mpz_cmp_ui(digit, 2));
  
  mpz_clear(digit);
  cf_free(x);
  cf_free(conv);

  x = cf_new_const(sqrt2);
  cf_t mob;
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
