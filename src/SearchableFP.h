#ifndef INCLUDED_UMDFP_FP_H
#define INCLUDED_UMDFP_FP_H 1
#define FINGERPRINT_SIZE 2048

#include <bitset>
#include <vector>
#include <fstream>
#include <cmath>

template<size_t N> static void saveBitSet(std::ostream& os, const std::bitset<N>& bv)
{
    unsigned char currentByte = 0;
    for (int i = 0; i < N; i++)
    {
        if (bv.test(i)) currentByte |= (1 << (i % 8)); // Set the bit in the current byte
        
        // Every 8th bit, or at the very end, write the byte
        if ((i + 1) % 8 == 0 || (i + 1) == N) {
            os.put(static_cast<char>(currentByte));
            currentByte = 0;
        }
    }
}

template<size_t N> static std::bitset<N> loadBitSet(std::istream& ist)
{
    std::bitset<N> ret;
    char read_byte;
    unsigned char bits;
    for(int i=0;i<N;)
    {
        ist >> read_byte;
        bits=static_cast<unsigned char>(read_byte);
        for(int j=0;j<8;j++) ret[i++]=(bits & (1 << j));
    }
    return ret;
}

template<size_t N> float tanimoto(std::bitset<N> b1, std::bitset<N> b2)
{
    std::bitset<N> un=(b1 | b2), isc=(b1 & b2);
    return ((float)isc.count())/un.count(); 
}

template<size_t N> float cosine(std::bitset<N> b1, std::bitset<N> b2)
{
    std::bitset<N> isc=(b1 & b2); // A.B
    float den=::sqrt((float)(b1.count()*b2.count()));
    return (float)(isc.count())/den;
}


#endif
