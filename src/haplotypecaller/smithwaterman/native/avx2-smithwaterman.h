#include <stdint.h>
#include "smithwaterman_common.h"

#undef SIMD_ENGINE
#define SIMD_ENGINE avx2

#include "avx2-functions.h"
#include "PairWiseSW.h"