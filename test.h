#define EXPECT(condition) \
  if (condition); else fprintf(stderr, "%s:%d: FAIL\n", __FILE__, __LINE__)
