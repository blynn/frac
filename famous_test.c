// Test famous continued fraction expansions are correct.
// as well as expansions based on Taylor series.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "cf.h"
#include "test.h"

int main() {
  CF_NEW_EXPECT_DEC(cf_new_e, "2718281828");
  CF_NEW_EXPECT_DEC(cf_new_pi, "31415926535897932384");
  CF_NEW_EXPECT_DEC(cf_new_tan1, "15574077246549022305");

  CF_NEW_EXPECT_DEC(cf_new_sin1, "08414709848078965066");
  CF_NEW_EXPECT_DEC(cf_new_cos1, "05403023058681397174");
  return 0;
}
