# Analyzing the Coronavirus (COVID-19) Genome

The coronavirus, also known as the Wuhan Coronavirus, `COVID-2019`  or `nCoV2019`, is a highly infectious disease that has to date
infected nearly 90,000 people around the world. Multiple nations across the globe have imposed major lockdowns and quarantines
to help contain it's spread. `nCoV2019` is believed to have originated in Wuhan, China, but has since spread to dozens
of other countries. South Korea and Italy have spawned new epicenters of coronavirus outbreaks, with thousands of cases
each. 

As the amount of cases and fatalities continues to rise, the question remains, what exactly is the coronavirus? There are
many different types of coronavirus. Some are close to harmless, while some are even more lethal than `nCoV2019`. But, that doesn't
tell us much about what makes this virus *different*. It has infected over 10 times the amount of people that SARS did, and has killed 
several times as many people. Yet no one knows what gives this virus an extra powerful infectious ability.

What we do know though, is the sequence of it's genome. The `nCoV2019` is an RNA virus. However, researchers from the 
Chinese Center for Disease Control and Prevention, Beijing, China [published a sequenced, reverse transcribed DNA genome](https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(20)30251-8/fulltext) of the
coronavirus. The `nCoV2019`'s genome is around 8,000 base pairs long. Far, far shorter than complex organisms like a human or rat. Yet, still powerful
enough to infect thousands of people across the entire world.

Over the course of this text, we will discuss the application of high performance search methods on DNA sequences. The abstraction
of DNA in the context of a programming language will also be discussed. This article will mainly focus on using `C` and `C++`, along
with `python` to analyze and search long sequences of DNA base pairs.

## DNA

Before, we start looking at the coronavirus genome harvested from patients in Wuhan, China, let's first understand
what DNA is. Deoxyribonucleicacid (DNA) is double stranded ladder-like structure composed of 4 different types of bases. They are
Thymine, Adenine, Guanine, and Cytosine. In sequencing, they are more commonly written as `T`, `A`, `G`, and `C` respectively. A DNA sequence
is represented as a finite sequences of bases, listed in the order in which they appear in the actual DNA. In a programming language, a
DNA sequence can be represented as a string only containing four possible characters, `'T', 'A', 'G', 'C'`.

In theory, the following can represent a string of DNA:

```c
const char* DNA = "ATGGCATGC";
```

There are a few issues with using a `char` type to represent a DNA base. First, a `char` can be any value between `-128` and 
`127`. Negative values definitvely aren't needed for genomic analysis, and we only need four values to represent the four possible
bases. Second, using a `const char*` literal forces each base to be represented by it's `ASCII` value. This value is not
in order, such as `'A'` being 65, yet `'G'` being 71. This prevents the bases from being uses as *indexes*, as they are not adjacent
to one antoher.

Lastly, it's more effecient to deal with data that's the size of a machine word, rather than data smaller than that.
For example, on most C/C++ compilers, fields in a struct with a size smaller than a machine word will be padded:

```c

struct foo {
   int a;
   char c;
};

printf("%zu\n", sizeof(struct foo));
```

The above example normally prints `8`, because it's more effecient to have the field sizes of `struct foo` align
with by 32bits than 8bits and be a total size of `5`.

Lastly, C strings are terminated by a null character `'\0'`. For DNA, that would be invalid, 
as `0` is a valid DNA sequence member. A single integer cannot terminate a DNA sequence,

### Types

*The code in this article is available at [this repository](https://github.com/jweinst1/ncov2019-analyzer)*

In order to overcome the shortcomings of the `char` type, an enumeration type wrapped in a `struct` 
may be used, with a case representing each base:

```cpp
struct DNA {
	enum Base {
		A = 0,
		C = 1,
		G = 2,
		T = 3
	};
};
```

Using an `enum` to represent DNA allows DNA data to be handled with type safety, as well as
having the other benefits described previously. Simiarly to C-strings, one can handle sequences of DNA
via pointers, such as `DNA::Base*` , or `const DNA::Base*`.

Given this type, a dna sequence could be constructed by doing

```cpp
DNA::Base* dna = new DNA::Base[3]{DNA::A, DNA::G, DNA::C};
// usage ... //
delete[] dna;
```

However, in order to effeciently iterate and search sequences of DNA, we need the proper way
to convert text data to dna, and data structures to ensapculate `DNA::Base*`.

### Conversion

Next, we need a way to convert dna textual data into the
`DNA::Base` type. Assuming the dna data is utf-8 or ASCII encoded, a function such as the following can work:

```cpp
struct DNA {
	enum Base {
		A = 0,
		C = 1,
		G = 2,
		T = 3
	};
    
    static void fromCStr(Base* dest, size_t dest_size, const char* src);
}

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
```

This approach assumes that dna data may either be upper case or lower case. Most `C++` compilers will
generate a jump table from a switch statement that has 5 or more cases. Therefore, the complexity of converting
`const char*` data into `const DNA::Base*` should only be linear with respective to the length of the
text data. This newly constructed sequence of DNA could be iterated and searched similar to how
strings can be searched, except with only four possible values at any index, `A, G, C, T`. To make
this more straight forward, we need data structures that are able of encapsulating dna. 

### Containers

The simplest dna data structure that can be used to hold dna is a *slice*. A slice is an object
which has a static size, and owns a chunk of memory with `Base::DNA` written to it. A slice can be
represented as:

```cpp
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
	
	const DNA::Base* dna() const 
	{
        return _dna;
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

```

In the above approach, either an *empty* `DNASlice` can be constructed, or a slice from a
C-string. An empty slice may be useful to represent a lack of matches from a search being run on DNA. The
`DNA::fromCStr()` function is reused here for the `const char*` constructor. Even though dna data
can now be encapsulated, we cannot visualize it. To do that, a customized printing method is needed to work on
`DNA::Base*` memory.

```cpp
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
``` 

This method directly prints characters to `stdout`. Since the `DNA::Base` enum has only 4 cases,
we do not need any `default` statement. One alternative would be to construct a `const char*` that *represents*
the base pairs of the DNA. This method uses less memory and doens't require calls to `new` and `delete`.

Next, DNA needs to be searched for some desired subsequence of dna. This can be accomplished similarly
to how one would search a string. Normally, this could just be done by calling every address in
the string and the target subsequence with `strcmp(str1, str2) == 0`. However, since the goal here is to illustrate the use of 
high performance methods, we can use a basic finite state machine to achieve better runtime performance.


```cpp
bool DNA::contains(const DNA::Base* dna1, 
                   size_t size1, 
                   const DNA::Base* dna2, 
                   size_t size2)
{
    bool mode = false;
    const DNA::Base* end1 = dna1 + size1;
    const DNA::Base* end2 = dna2 + size2;
    const DNA::Base* matcher = dna2;
    if (size2 > size1)
        return false;
    while (dna1 != end1) {
        if (mode) {
            if (matcher == end2)
                return true;
            else if (*matcher == *dna1)
                ++matcher;
            else {
                matcher = dna2;
                mode = false;
            }
        } else {
            if (*matcher == *dna1) {
                mode = true;
                ++matcher;
            }
        }
        ++dna1;
    }
    // Check once here incase it's at the end
    if (matcher == end2)
        return true;
    return false;
}
```

The above implementation works by detecting for a match of `dna2` at any starting position in `dna1`. 
It aliases `dna2` through `matcher`, this allows `matcher` to be reset in the case that only a partial
match is found. Before any iteration we know that if `size2 > size1` is true, it's impossible for
`dna2` to be contained in `dna1`.

The boolean variable: 

```cpp
bool mode = false;
```

is used to determine when the function is in a *state* of matching and when it's in an
*initial* state, still looking for the first base to match. If in the initial state, we only
transition into the matching state if we find the first base. Now, in the matching state, the following
conditons apply:

```cpp
    if (matcher == end2)
        return true;
    else if (*matcher == *dna1)
        ++matcher;
    else {
        matcher = dna2;
        mode = false;
    }
```

Either we find match the last base of `dna2`, and thus find a *full* match, we proceed further
toward a full match, or we fail to match a base before reaching the end, causing a reset
out of the matching state.

### Search

Now that conversion, iteration, and encapsulation of DNA have been discussed and demonstrated, we can showcase
how to search DNA. Previously, using an FSA was shown as a method of linear search for a particular subsequence of 
DNA. That approach is useful, since it can count multiple occurences of some subsequence of dna in a larger
sequence. However, it lacks the performance to be able to scan and identify against many, many sequences of DNA.

An even mor effecient approach is to use a trie structure that's specialized for dna. Normally,
although a trie has very fast lookup time, it's downside is using a lot of memory due to the
amount and size of each node that is needed. Typically, a trie node could appear like this.

```cpp

struct TrieNode {
    TrieNode* childNodes[128];
    // ... other attributes //
};
```

Here, the node has 128 child node slots because the maximum value for the `char` type in
C is *usually* 127. However, that might differ depending on the value of `SCHAR_MAX` defined in `<limits.h>`. 
In terms of space complexity, a trie using a node defined like this will take up a lot of memory with a
decent amount of keys inserted. In general, the space and memory used by a trie correlatr with the
range of characters it supports in it's keys.

#### DNA Tries

In the case of DNA, there are only 4 possible characters! `A, G, C, T`. Therefore, a dna trie node
would look like this.

```cpp
struct DNANode {
    DNANode* childNodes[4];
    // other properties ... //
};
```

For any given base at any given position in a DNA sequence, there are only four possible values
for the next base. Thus, this principal can be used to make a `DNA::Base` function as the index 
in the child nodes of a `DNANode`. The child nodes can then be accessed like

```cpp
DNA::Base* dna = new DNA::Base[3]{ DNA::A, DNA::A, DNA::G}
DNANode node;
node.childNodes[dna[1]]; // Use dna to get child node.
delete[] dna;
```

If a node only requires around 32-40 bytes, this makes it far more space effecient than
a node covering the entire signed range of `char`. Next, we can choose what attirbutes and values our trie will map
with it's keys. There are three important types of data when searching dna. One is existence, such as
if a subsequence of dna is contained in a larger sequence. Another is counting the number of times
a certain subsequence appears. Lastly, dna can be searched for the areas that a subsequence appears the most.

For the purposes of this article, let's use a trie which keeps a count of the DNA sequences that
match that particular node. Assume that we may only want to monitor for particular sub regions of DNA,
and thus small prefixes of DNA such as `AG` or `CT` are not meaningful. The trie will then have a count of `-1`
for any node initially, and upon insertion, will increment to 0, and so on. This allows us to track only
important sequences in the trie.

The core trie class will look like this:

```cpp
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
```

First, the node class should have a flexible constructor that can take 0 to 4 arguments
as child nodes. The node should also possess a method to determine if it's a leaf or not. The core functionality
of the node class comes from the `insert` and `find` methods. A `remove` method can also be implemented,
but it's not critical for the scope of this text.

To insert into a trie, we need to recursively traverse the trie, and create new nodes for whatever base 
paths do not yet exist in the trie. Once the end of the input dna is reached, the count of the node is incremented by 1. 
There is the option to potentially view the input dna as a nested sequence, and count all subsequences contained within it.
Such that if the input is `ACCG`, we would increment `A`, `AC`, `ACC`, not just `ACCG`.

```cpp
void DNANode::insert(const DNA::Base* dna, size_t size, bool nested)
{
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
```

In the recursion scheme, the base case is `size == 0`, because this indicates the end of the dna input. `++count`
is used instead of `count++` for a slight performance boost.

Next, the trie must be able to lookup the counts associated with each DNA sequence inserted into the trie.
The approach is similar to `insert`, but returns a `long`:

```cpp
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
```

If `find` returns `-1`, that means the dna sequence does not exist in the trie. If it returns `0`,
it means the sequence does exist in the trie, but has not been inserted more than once. Subsequent return values
indicate the count of the sequence tracked within the trie. Although the `TrieNode` can be used directly, 
it's more effecient to wrap it in another class that manages a root node where all incoming insertions and lookups 
are routed to. Such a class could look like this:

```cpp

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

```

In the above class, `_count` keeps track of the number of sequences, yet not unique sequences, inserted
into the trie. Otherwise, there would be no fast way to count the number of strings counted. The `<<` operator
is used for insertion to make the class feel similar to stream class. The intention is for the trie to be able to
insert a large number of DNA chunks.

The DNA Trie is best suited for when the important sequences one is looking for are known. That way, varying 
chunks of DNA can go in, and `find()` can be used on particular, important sequences.

## The COVID-19 Genome

Finally, we have discussed the tools and techniques to analyze DNA in a high performance context.
Next, let's talk about the actual coronavirus genome.

This coronavirus, like most, are RNA viruses. Yet, DNA can be transcribed from RNA, to make analysis more consistent.
The genome of `COVID-19` is about 29kb on disk. This is not representative of any COVID-19 in existence, 
this particular sequence was taken from patients in Wuhan, China. Viruses can evolve as infections spread
overtime.

*The file containing the genome can be found [here](https://github.com/jweinst1/ncov2019-analyzer/blob/master/data/COVID-19-genome.txt)*

*Credits for the sequence are given to  Credits to Wu,F., Zhao,S., Yu,B., Chen,Y.M., Wang,W., Song,Z.G., Hu,Y., ao,Z.W., Tian,J.H., Pei,Y.Y., Yuan,M.L., Zhang,Y.L., Dai,F.H., Liu,Y., Wang,Q.M., Zheng,J.J., Xu,L., Holmes,E.C. and Zhang,Y.Z.*

The first step in analyzing the genome is deciding what kind of patterns are important to look for, and how big are those patterns.
Say for instance, we want to see how many times the pattern `ACGGTTCCAAT` occurs. We could read the COVID-19 genome and slice it
every 11 bases, and feed it into an instance of a DNA Trie.

### Formulating a search

Based on the techniques discussed previously, we can now form a search of what we want to look for in the COVID-19
genome. Ideally, it's desirable to search for one or more sequences at a time. A solution is needed
that will accept an arbitrary amount of `const char*` strings, convert them to `DNA::Base*`, and insert them
into a `DNATrie`. For simplicity, the `main()` function and it's variable arguments can be used. We will
also need a class that can be used within a `vector<>` to hold potentially many search sequences and store the
results of each search, such as:

```cpp
struct GenomeArgument {
    GenomeArgument(const char* genseq): seq(genseq),
                                        sSize(std::strlen(seq)),
                                        fSize(0)
                                        {}
    void check(DNATrie& trie)
    {
        fSize = trie[seq];
    }
    
    static void populate(char const** args, int size, std::vector<GenomeArgument>& gens);
    
    const char* seq;
    size_t sSize;
    long fSize;
};

void GenomeArgument::populate(char const** args, int size, std::vector<GenomeArgument>& gens)
{
    while (size--) {
        GenomeArgument garg(*args++);
        gens.push_back(std::move(garg));
    }
}
```


