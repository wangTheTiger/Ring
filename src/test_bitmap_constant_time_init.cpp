#include<iostream>
#include <sdsl/init_array.hpp>

typedef uint16_t word_t;

int main(int argc, char **argv)
{
    int max_size = 300;
    sdsl::initializable_array<word_t> D_array(4*(max_size+1), 0);
    cout << " D_array initialized in constant time";
}