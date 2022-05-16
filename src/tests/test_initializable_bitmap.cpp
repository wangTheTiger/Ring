
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
    //std::cout << bitmap.atPos(5) << std::endl;
    bitmap[5] |= 1ULL << (5 % 64);
    //std::cout << bitmap.atPos(5) << std::endl;
    //std::cout << bitmap.atPos(0) << std::endl;
    //std::cout << bitmap.atPos(1) << std::endl;
    bitmap[128] |= 1ULL << (128 % 64);
    bitmap[257] |= 1ULL << (257 % 64);
    bitmap[511] |= 1ULL << (511 % 64);
    bitmap[512] |= 1ULL << (512 % 64);
    /*
    Then, we print all the bits out :-).
    */
    for (int i = 0 ; i < MAX; i++){
        std::cout << bitmap.atPos(MAX - i);
    }
    std::cout << " " << std::endl;
}