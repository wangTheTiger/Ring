
#include <sdsl/init_array.hpp>
#include <sdsl/int_vector.hpp>
#include<iostream>
int main()
{
    uint64_t MAX = 10;
    //>>Bit Vector testing
    sdsl::bit_vector bv(MAX, 0);
    for (uint64_t i = 0; i < MAX; i++)
    {
        bv[i] = 1;
    }
    //<Bit Vector testing
    //Building an init array of bit_vectors. (Not what I want but just to demonstrate an idea).
    sdsl::initializable_array<sdsl::bit_vector> D_array(MAX, sdsl::bit_vector(MAX,0));
    D_array[0] = sdsl::bit_vector(MAX,0);
    D_array[1] = sdsl::bit_vector(MAX,1);
    D_array[2] = sdsl::bit_vector(MAX,1);
    std::cout << "D_array " << ((float)D_array.size_in_bytes()) << " bytes." << std::endl;
    std::cout << "D_array.atPos(0) = " << D_array.atPos(0) << std::endl;
    std::cout << "D_array.atPos(1) = " << D_array.atPos(1) << std::endl;
    std::cout << "D_array.atPos(2) = " << D_array.atPos(2) << std::endl;
    //The basic data structure at the hardware level of mainstream CPUs is a byte.
    // Necesito algo asi, un arreglo inializable que reciba un tipo que tenga un solo bit. sdsl::initializable_array<1> D_array2(MAX, 1); <-- no se puede por comentario de arriba.
    sdsl::initializable_array<uint8_t> D_array2(MAX, 0);
    std::cout << "Initializable array of uint8_t (8 bytes) " << ((float)D_array2.size_in_bytes()) << " bytes = 8 bytes per cell x 10 cells + 8 bytes pointer." << std::endl;
}