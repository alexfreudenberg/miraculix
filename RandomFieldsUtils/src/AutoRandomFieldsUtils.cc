#include "AutoRandomFieldsUtils.h"

const char

*LA_NAMES[LA_LAST + 1] = {"auto", "intern", "R", "GPU"},

*PIVOT_NAMES[PIVOT_LAST + 1] = {"none", "auto", "do", "idx", "undefined"},

*INSTALL_NAMES[INSTALL_LAST + 1] = {"none", "ask", "install", 
				    "sse", "sse2", "sse3", "ssse3", "avx", "avx2",
				    "gpu"};
