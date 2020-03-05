#ifndef EXAMPLES_TEST_OBJECT_H
#define EXAMPLES_TEST_OBJECT_H

#include <cstdlib>
#include <cstdio>
#include <cstring>

class TestObject {
public:
    void clear()
    {
        _failures = 0;
    }
    
    template<class T>
    void eq(const T& lfs, const T& rfs)
    {
        _failures += (unsigned)(!(lfs == rfs));
    }
    
    void isTrue(bool expres)
    {
        _failures += !expres ? 1 : 0;
    }
    
    void isFalse(bool expres)
    {
        _failures += expres ? 1 : 0;
    }
    
    template<class T>
    void lt(const T& lfs, const T& rfs)
    {
        _failures += (unsigned)(!(lfs < rfs));
    }

    template<class T>
    void gt(const T& lfs, const T& rfs)
    {
        _failures += (unsigned)(!(lfs > rfs));
    }
    
    template<class T>
    void notNull(const T* obj)
    {
        _failures += (unsigned)(obj == nullptr);
    }
    
    void finish()
    {
        std::printf("Got %u failures\n", _failures);
        std::exit(_failures ? 3 : 0);
    }
private:

    
    unsigned _failures;
};

#endif // EXAMPLES_TEST_OBJECT_H