#define EXPECT(condition) \
  if (condition); else fprintf(stderr, "%s:%d: FAIL\n", __FILE__, __LINE__)

#define CF_EXPECT_DEC(cf, str) \
  cf_expect_dec(cf, str, __FILE__, __LINE__)

#define CF_NEW_EXPECT_DEC(cf, str) \
  cf_new_expect_dec(cf, str, __FILE__, __LINE__)

void cf_expect_dec(cf_t x, char *result, char *filename, int line);

void cf_new_expect_dec(cf_t (*cf_new_fn)(), char *result,
    char *filename, int line);
