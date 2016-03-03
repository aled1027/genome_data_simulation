
#include <stdint.h>
#include <stdio.h>
#include <inttypes.h>
#include <assert.h>
#include <math.h>

void print64(uint64_t x)
{
    printf("%" PRIu64, x);
}

uint64_t hash(uint64_t key) {
    // https://naml.us/blog/2012/03/inverse-of-a-hash-function
    key = (~key) + (key << 21); // key = (key << 21) - key - 1;
    key = key ^ (key >> 24);
    key = (key + (key << 3)) + (key << 8); // key * 265
    key = key ^ (key >> 14);
    key = (key + (key << 2)) + (key << 4); // key * 21
    key = key ^ (key >> 28);
    key = key + (key << 31);
    return key;
}

int main(int argc, char *argv[]) 
{
    assert(argc && argv);
    uint64_t start = 0;
    uint64_t end = pow(2,10);
    for (uint64_t i = start; i < end; ++i) {
        print64(i);
        printf(", ");
        print64(hash(i));
        printf("\n");
    }
    return 0;
}
