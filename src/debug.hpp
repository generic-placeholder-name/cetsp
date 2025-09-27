#pragma once

#include <iostream>

// Compile with -DDEBUG to enable debug logging
#if DEBUG
    #define DBG(x) do { std::cerr << x << std::endl; } while (0)
    #define DBG_NOENDL(x) do { std::cerr << x; } while (0)
#else
    #define DBG(x) do {} while (0)
    #define DBG_NOENDL(x) do {} while (0)
#endif
