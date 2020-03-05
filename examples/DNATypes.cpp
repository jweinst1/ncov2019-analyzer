#include <algorithm>
#include <cstdio>
#include <cstring>
#include "TestObject.h"

static TestObject gTester;

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
    static void print(const Base* dna, size_t dna_size);
    static bool contains(const Base* dna1, 
                         size_t size1, 
                         const Base* dna2, 
                         size_t size2);
};

void DNA::fromCStr(DNA::Base* dest, size_t dest_size, const char* src)
{
    while (dest_size-- && *src) {
        switch (*src) {
            case 'a':
            case 'A':
                *dest++ = A;
                break;
            case 'c':
            case 'C':
                *dest++ = C;
                break;
            case 'g':
            case 'G':
                *dest++ = G;
                break;
            case 't':
            case 'T':
                *dest++ = T;
                break;
            default:
                // Ignore non DNA chars
                break;
        }
        ++src;
    }
}

void DNA::print(const DNA::Base* dna, size_t dna_size)
{
    while (dna_size--) {
        switch (*dna++) {
             case DNA::A:
                 std::putc('A', stdout);
                 break;
             case DNA::C:
                 std::putc('C', stdout);
                 break;
             case DNA::G:
                 std::putc('G', stdout);
                 break;
             case DNA::T:
                 std::putc('T', stdout);
                 break;
        }
    }
}

bool DNA::contains(const DNA::Base* dna1, 
                   size_t size1, 
                   const DNA::Base* dna2, 
                   size_t size2)
{
    int mode = 0;
    const DNA::Base* end1 = dna1 + size1;
    const DNA::Base* end2 = dna2 + size2;
    const DNA::Base* matcher = dna2;
    if (size2 > size1)
        return false;
    while (dna1 != end1) {
        if (mode) {
            
        } else {
            if (*matcher == *dna1) {
                
            }
        }
        ++dna1;
    }
    return false;
}

static void testfromCStr()
{
    size_t data_size = 10;
    DNA::Base* data = new DNA::Base[data_size];
    DNA::fromCStr(data, data_size, "AATG");
    std::printf("Got the DNA ");
    DNA::print(data, data_size);
    std::puts(" .");
    gTester.eq(data[0], DNA::A);
    gTester.eq(data[1], DNA::A);
    gTester.eq(data[2], DNA::T);
    gTester.eq(data[3], DNA::G);
    delete[] data;
}

class DNASlice {
public:
    DNASlice();
    explicit DNASlice(const char* data);
    ~DNASlice();
    
    bool empty() const
    {
        return !_size;
    }
	
	size_t size() const
	{
		return _size;
	}
	
	bool operator==(const DNASlice& other) const;
	bool operator!=(const DNASlice& other) const
	{
		return !(*this == other);
	}
private:
    size_t _size;
    DNA::Base* _dna;
};

DNASlice::DNASlice(): _size(0), _dna(nullptr) {}

DNASlice::DNASlice(const char* data): _size(std::strlen(data)),
                                      _dna(new DNA::Base[_size])
{
    DNA::fromCStr(_dna, _size, data);
}

DNASlice::~DNASlice()
{
    delete[] _dna;
}

bool DNASlice::operator==(const DNASlice& other) const
{
	if (_size != other._size)
		return false;
	return std::equal(_dna, _dna + _size, other._dna, other._dna + _size);
}

class DNAView {
public:
   DNAView(const DNA::Base* dna, size_t size): _size(size),
                                               _dna(dna)
                                               {}
    size_t size() const
    {
        return _size;
    }
    
private:
   size_t _size;
   const DNA::Base* _dna;   
};

static void testDNASlice()
{
	static const char* TESTDNA = "AGCCT";
	DNASlice slc(TESTDNA);
	DNASlice slc2("ACGTT");
	DNASlice foo("ACGTT");
	gTester.eq(foo, slc2);
	gTester.isFalse(slc == slc2);
	
}

int main(int argc, char const* argv[]) {
    gTester.clear();
    testfromCStr();
	testDNASlice();
    gTester.finish();
	return 0;
}