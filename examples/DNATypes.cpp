#include <algorithm>
#include <exception>
#include <cstdio>

/**
 * Acts as a scalar type to represent DNA bases,
 * ensure the cases are sequential 
 */
struct DNA {
	enum Base {
		A = 0,
		C = 1,
		G = 2,
		T = 3
	};
    
    static void fromCStr(Base* dest, size_t dest_size, const char* src);
};

void DNA::fromCStr(DNA::Base*, size_t dest_size, const char* src)
{
    while (dest_size && *src) {
        switch (*src) {
            case 'a':
            case 'A':
                break;
            case 'c':
            case 'C':
                break;
            case 'g':
            case 'G':
                break;
            case 't':
            case 'T':
                break;
            default:
                break;
        }
    }
}

int main(int argc, char const* argv[]) {
	return 0;
}