#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>

/**
 * Genome analyzer for COVID19
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
    
    DNAChunk(DNA::Base* data, size_t size): dsize(size),
                                           dna(new DNA::Base[dsize])
   {
       std::copy(data, data + size, dna);
   }                                       
    ~DNAChunk()
    {
        delete[] dna;
        dsize = 0;
    }
    
    DNAChunk(DNAChunk&& other): dsize(other.dsize),
                                dna(other.dna)
                                {}
    
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
                
    ~DNANode();
				
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
    
    void insert(const DNA::Base* dna, size_t size, bool nested = false);
    long find(const DNA::Base* dna, size_t size);
    void remove(const DNA::Base* dna, size_t size);
	
	long count;
    DNANode* childNodes[4];
};

static char _dnaBaseToChar(DNA::Base b)
{
    switch(b) {
        case DNA::A: return 'A';
        case DNA::C: return 'C';
        case DNA::G: return 'G';
        case DNA::T: return 'T';
    }
    return '\0';
}

void DNANode::insert(const DNA::Base* dna, size_t size, bool nested)
{
#ifdef TRIE_TEST_DEBUG
        std::printf("Inserting at base: %c, size %zu\n", _dnaBaseToChar(*dna), size);
#endif
    if (!size) {
        ++count;
        return;
    }
    DNANode** nextNode = &childNodes[*dna];
    if (*nextNode == nullptr) {
        *nextNode = new DNANode();
    }
    if (nested)
        ++count;
    (*nextNode)->insert(dna + 1, size - 1);
}

long DNANode::find(const DNA::Base* dna, size_t size)
{
    if (!size)
        return count;
    DNANode* nextNode = childNodes[*dna];
    if (nextNode == nullptr)
        return -1;
    else
        return nextNode->find(dna + 1, size - 1);
}

void DNANode::remove(const DNA::Base* dna, size_t size)
{
    if (!size) {
        count = -1;
        return;
    }
    DNANode* nextNode = childNodes[*dna];
    if (nextNode != nullptr)
        nextNode->remove(dna + 1, size - 1);
}

DNANode::~DNANode()
{
    if(childNodes[0] != nullptr) delete childNodes[0];
    if(childNodes[1] != nullptr) delete childNodes[1];
    if(childNodes[2] != nullptr) delete childNodes[2];
    if(childNodes[3] != nullptr) delete childNodes[3];
    childNodes[0] = nullptr;
    childNodes[1] = nullptr;
    childNodes[2] = nullptr;
    childNodes[3] = nullptr;
}

class DNATrie {
public:
   DNATrie(): _count(0){}
   ~DNATrie(){}
   
   const DNANode& getRoot() const
   {
       return _root;
   }
   
   size_t getCount() const
   {
       return _count;
   }
   
   DNATrie& operator<<(const char* input)
   {
       DNAChunk chk(input);
       _root.insert(chk.dna, chk.dsize);
       ++_count;
       return *this;
   }
   
   DNATrie& operator<<(const DNAChunk& input)
   {
       _root.insert(input.dna, input.dsize);
       ++_count;
       return *this;
   }
   
   long operator[](const char* input)
   {
       DNAChunk chk(input);
       return _root.find(chk.dna, chk.dsize);
   }
   
   
private:
   size_t _count;
   DNANode _root;
};
 
static void usage_print()
{
    std::puts("____HELP_____");
    std::puts("Usage:");
    std::puts("$ covid19 <path_to_genome> <sequence_to_search>");
    std::puts("_____________");
}

static TrieNode gTrie;

int main(int agrc, char const* argv[])
{
    FILE* gfp;
    char* readBuf;
    size_t searchSize;
    if (argc != 3) {
        usage_print();
        std::exit(1);
    }
    // test make sure file can be opened
    gfp = std::fopen(argv[1], "r");
    if (gfp == NULL) {
        fprintf("ERR: Genome file at '%s' cannot be opened.\n", argv[1]);
        std::exit(2);
    }
    searchSize = std::strlen(argv[2]);
    readBuf = new char[searchSize + 1](0);
    while (fgets(readBuf, searchSize, gfp) != NULL) {
        // Process the read chunk of genome.
    }
    
    
    return 0;
}