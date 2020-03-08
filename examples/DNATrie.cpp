#include <algorithm>
#include <cstdio>
#include <cstring>
#include <cassert>
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

struct DNAChunk {
    DNAChunk(const char* string): dsize(strlen(string)),
                                  dna(new DNA::Base[dsize])
   {
       DNA::fromCStr(dna, dsize, string);
   }
    ~DNAChunk()
    {
        delete[] dna;
        dsize = 0;
    }
    
    size_t dsize;
    DNA::Base* dna;
};

struct DNANode {
				
	DNANode(DNANode* d1 = nullptr,
	        DNANode* d2 = nullptr,
			DNANode* d3 = nullptr,
			DNANode* d4 = nullptr): count(-1),
			childNodes{d1, d2, d3, d4}
				{}
				
	DNANode* operator[](const DNA::Base base)
	{
		return childNodes[base];
	}
	
	bool isLeaf() const
	{
		return childNodes[0] == nullptr &&
		       childNodes[1] == nullptr &&
			   childNodes[2] == nullptr &&
			   childNodes[3] == nullptr;
	}
    
    void insert(const DNA::Base* dna, size_t size);
    long find(const DNA::Base* dna, size_t size);
	
	long count;
    DNANode* childNodes[4];
};

void DNANode::insert(const DNA::Base* dna, size_t size)
{
    assert(size);
    if (size == 1) {
        count = count == -1 ? 0 : count + 1;
        return;
    }
    DNANode** nextNode = &childNodes[*dna];
    if (*nextNode == nullptr) {
        *nextNode = new DNANode();
    }
    (*nextNode)->insert(dna + 1, size - 1);
}

long DNANode::find(const DNA::Base* dna, size_t size)
{
    assert(size);
    if (size == 1)
        return count;
    DNANode* nextNode = childNodes[*dna];
    if (nextNode == nullptr)
        return -1;
    else
        return nextNode->find(dna + 1, size - 1);
}

int main(int argc, char const* argv[])
{
    gTester.clear();
    
    gTester.finish();
    return 0;
}