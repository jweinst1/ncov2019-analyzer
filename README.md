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
text data. 