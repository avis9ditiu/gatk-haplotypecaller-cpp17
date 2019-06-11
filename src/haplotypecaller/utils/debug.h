#ifndef DEBUG_H
#define DEBUG_H

#include <cstdio>
#ifdef _DEBUG
#define DBG(format, args...) printf("[%s:%d] "format"\n", __FILE__, __LINE__, ##args)
#else
#define DBG(args...)
#endif

#endif //DEBUG_H
