
/*! \file test_initializable_bitmap.hpp
    \brief A simple test created for init_bitmap.hpp (SDSL) which implements an initializable bitmap based on https://users.dcc.uchile.cl/~gnavarro/ps/acmcs14.pdf and init_array.hpp.
    Three functions initialize a B array cell, they are operators +=, -= and [](uint64_t _i).
    TODO: move it to SDSL :-)
    \author Fabrizio Barisione
*/
#include <sdsl/init_bitmap.hpp>
#include <sdsl/int_vector.hpp>
#include<iostream>
int main()
{
    uint64_t MAX = 512;
    sdsl::initializable_bitmap bitmap(MAX, 0);
    std::cout << "Initializable bitmap, size in bytes = "<< ((float)bitmap.size_in_bytes()) << std::endl;
    cout << "bitmap: ";
    /*
    Turing on bits at positions 5, 128, 257, 511 and 512.
    We get the memory address of the block as a whole (word size) and then we set as 1 the bit within the block at offset i & 64.
    */

    /*OLD WAY, without += (set) and -= (clear) bit operators.
    bitmap[5] |= 1ULL << (5 % 64);
    bitmap[128] |= 1ULL << (128 % 64);
    bitmap[257] |= 1ULL << (257 % 64);
    bitmap[511] |= 1ULL << (511 % 64);
    bitmap[512] |= 1ULL << (512 % 64);
    */

    //Setting bits: the new way.
    bitmap+=5;
    bitmap+=128;
    bitmap+=257;
    bitmap+=511;
    bitmap+=512;
    //clearing bits.
    bitmap-=512;
    bitmap-=5;
    bitmap-=0;
    /*
    Then, we print all the bits out :-).
    */
    for (int i = 0 ; i < MAX; i++){
        std::cout << bitmap.atPos(MAX - i);
    }
    std::cout << " " << std::endl;
}