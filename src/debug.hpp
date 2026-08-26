#pragma once

#include <iostream>

// Compile with -DDEBUG to enable debug logging
#if DEBUG
    #define DBG(x) do { std::cerr << x << std::endl; } while (0)
    #define DBG_NOENDL(x) do { std::cerr << x; } while (0)
#else
    // Keep expressions type-checked and variables visibly used while allowing
    // the optimizer to discard all disabled logging.
    #define DBG(x) do { if constexpr (false) { std::cerr << x << std::endl; } } while (0)
    #define DBG_NOENDL(x) do { if constexpr (false) { std::cerr << x; } } while (0)
#endif
