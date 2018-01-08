#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <algorithm>
#include <xmmintrin.h>
#include <pmmintrin.h>
using std::min;
using std::max;
typedef uint32_t DWORD;
#if 0
#define _aligned_malloc(a,b) aligned_alloc(b,a)
#else
#define _aligned_malloc(ptr, a,b) posix_memalign(ptr, b,a)
#endif
#define _aligned_free(a)