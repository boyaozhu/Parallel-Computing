#include <inttypes.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


#define ROTL32(x,y) ((x<<y)|(x>>(32-y)))
#define ROTR32(x,y) ((x>>y)|(x<<(32-y)))
#define ROTL24(x,y) ((x<<y)|(x>>(24-y)))
#define ROTR24(x,y) ((x>>y)|(x<<(24-y)))
#define ROTL16(x,y) ((x<<y)|(x>>(16-y)))
#define ROTR16(x,y) ((x>>y)|(x<<(16-y)))
#define ROTL8(x,y)  ((x<<y)|(x>>(8-y)))
#define ROTR8(x,y)  ((x>>y)|(x<<(8-y)))


void encrypt32( char * filename, uint32_t key ) {

    FILE * in, * out;
    char tempfilename[9999];
    uint32_t block;
    uint32_t a, b, c, d, magnitude, polarity, xor;

    srand(key);
    sprintf(tempfilename, "%s%s", filename, ".ecp");

    if ((in = fopen(filename, "rb")) == NULL) {

        fprintf(stderr, "ERROR: Could not open input file %s\n", filename);
        exit(1);
    }

    if ((out = fopen(tempfilename, "wb")) == NULL) {

        fprintf(stderr, "ERROR: Could not open output file %s\n", tempfilename);
        exit(1);
    }

    /* main loop for encoding - all 4 bytes are valid */
    while (     (a = fgetc(in)) != EOF
            &&  (b = fgetc(in)) != EOF
            &&  (c = fgetc(in)) != EOF
            &&  (d = fgetc(in)) != EOF
    ) {
        printf( "a = %c, b = %c, c = %c, d = %c\n", a, b, c, d );

        polarity = rand() % 2;
        magnitude = rand() % 32;

        block = (d << 24) | (c << 16) | (b << 8) | a;

        xor =   (rand() % 256 << 24) | (rand() % 256 << 16)
              | (rand() % 256 << 8)  | (rand() % 256);

        block ^= xor;

        if (polarity)
            block = ROTL32(block, magnitude);
        else
            block = ROTR32(block, magnitude);

        putc( block,       out );
        putc( block >> 8,  out );
        putc( block >> 16, out );
        putc( block >> 24, out );

        printf( "p = %d, mag = %d, xor = %d\n", polarity, magnitude, xor );
        printf( "da = %x, db = %x, dc = %x, dd = %x\n\n",
                block, block >> 8, block >> 16, block >> 24 );
    }

    /* cases to handle the end of file */
    if ( a != EOF && b != EOF && c != EOF && d == EOF ) {

        polarity = rand() % 2;
        magnitude = rand() % 24;

        block = (c << 16) | (b << 8) | a;

        block ^= (rand() % 256 << 16) | (rand() % 256 << 8) | (rand() % 256);

        if (polarity)
            block = ROTL24(block, magnitude);
        else
            block = ROTR24(block, magnitude);

        putc( block,              out );
        putc( block = block >> 8, out );
        putc( block = block >> 8, out );

    } else if ( a != EOF && b != EOF && c == EOF ) {

        polarity = rand() % 2;
        magnitude = rand() % 16;

        block = (b << 8) | a;

        block ^= (rand() % 256 << 8) | (rand() % 256);

        if (polarity)
            block = ROTL16(block, magnitude);
        else
            block = ROTR16(block, magnitude);

        putc( block,              out );
        putc( block = block >> 8, out );

    } else if ( a != EOF && b == EOF ) {

        polarity = rand() % 2;
        magnitude = rand() % 8;

        block = a;
        block ^= (rand() % 256);

        if (polarity)
            block = ROTL8(block, magnitude);
        else
            block = ROTR8(block, magnitude);

        putc(block, out);
    }

    fclose(in);
    fclose(out);
}


int main( int argc, char *argv[] ) {

    printf("x0r 32-bit File Encryption\n");

    if (argc < 3) {
        fprintf( stderr, "ERROR: No file(s) or encryption key supplied.\n" );
        fprintf( stderr, "USAGE: This program requires a filename " );
        fprintf( stderr, "and a key in the range [0,%u] as arguments!\n",
                 (unsigned) pow(2, sizeof(uint32_t) * 8) - 1);
        exit(1);
    }

    printf( "encryting file %s using key %s...\n", argv[1], argv[2] );
    printf( "To quit, press ctrl + c\n\n" );
    printf( "Status: \n" );

    encrypt32( argv[1], atol(argv[2]) );

    printf( "done!\n" );

    return 0;
}
