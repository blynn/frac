// Test famous continued fraction expansions are correct.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "cf.h"
#include "test.h"

int main() {
  mpz_t z;
  mpz_init(z);
  cf_t e, conv;
  e = cf_new_e();

  conv = cf_new_cf_to_decimal(e);
  char s[128];
  for (int i = 0; i <= 10; i++) {
    cf_get(z, conv);
    s[i] = mpz_get_ui(z) + '0';
  }
  s[10] = 0;
  EXPECT(!strcmp(s, "2718281828"));
  cf_free(conv);
  cf_free(e);

  cf_t pi;
  pi = cf_new_pi();
  conv = cf_new_cf_to_decimal(pi);
  for (int i = 0; i <= 20; i++) {
    cf_get(z, conv);
    s[i] = mpz_get_ui(z) + '0';
  }
  s[20] = 0;
  EXPECT(!strcmp(s, "31415926535897932384"));
  cf_free(conv);
  cf_free(pi);

  mpz_clear(z);
  return 0;
}
