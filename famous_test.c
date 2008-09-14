// Test famous continued fraction expansions are correct.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "cf.h"
#include "test.h"

void check_match(cf_t (*cf_new_fn)(), char *result) {
  mpz_t z;
  mpz_init(z);
  cf_t x, conv;
  x = cf_new_fn();
  conv = cf_new_cf_to_decimal(x);
  int len = strlen(result);
  char s[len + 1];
  for (int i = 0; i <= len; i++) {
    cf_get(z, conv);
    s[i] = mpz_get_ui(z) + '0';
  }
  s[len] = 0;
  EXPECT(!strcmp(s, result));

  cf_free(conv);
  cf_free(x);
  mpz_clear(z);
}

int main() {
  check_match(cf_new_e, "2718281828");
  check_match(cf_new_pi, "31415926535897932384");
  check_match(cf_new_tan1, "15574077246549022305");

  mpz_t z;
  mpz_init(z);
  cf_t e = cf_new_e();
  cf_t pi = cf_new_pi();
  mpz_t addarray[8];
  for (int i = 0; i < 8; i++) {
    mpz_init(addarray[i]);
    mpz_set_ui(addarray[i], 0);
  }
  mpz_set_ui(addarray[1], 1);
  mpz_set_ui(addarray[2], 1);
  mpz_set_ui(addarray[7], 1);
  cf_t b = cf_new_bihom(e, pi, addarray);
  cf_t conv = cf_new_cf_to_decimal(b);
  cf_get(z, conv);
  gmp_printf("e + pi = %Zd.", z);
  for (int i = 0; i <= 20; i++) {
    cf_get(z, conv);
    gmp_printf("%Zd", z);
  }
  printf("\n");

  mpz_clear(z);
  return 0;
}
