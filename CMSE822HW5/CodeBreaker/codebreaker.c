#include <inttypes.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#define ROTL32(x,y) ((x<<y)|(x>>(32-y)))
#define ROTR32(x,y) ((x>>y)|(x<<(32-y)))
#define ROTL24(x,y) ((x<<y)|(x>>(24-y)))
#define ROTR24(x,y) ((x>>y)|(x<<(24-y)))
#define ROTL16(x,y) ((x<<y)|(x>>(16-y)))
#define ROTR16(x,y) ((x>>y)|(x<<(16-y)))
#define ROTL8(x,y)  ((x<<y)|(x>>(8-y)))
#define ROTR8(x,y)  ((x>>y)|(x<<(8-y)))

int num_entries = 8;
const char * dictionary [] = {  "bird",
                                "campus",
                                "class",
                                "of",
                                "spring",
                                "sun",
                                "the",
                                "tree" };

// checks if a significant portion of the words resulting from
// the decryption with the tried key are found in the dictionary
int isValid( char *decoded, int len ) {

    int nmatches = 0;
    char * word;

    word = strtok( decoded, " ,.;-()\n\r" );

    while ( word != NULL ) {

        int flag = 0;

        for ( int i = 0; i < num_entries; ++i ) {
            if( !strcmp( word, dictionary[i] ) ) {
                flag = 1;
                break;
            }
        }

        if (flag) {  nmatches += strlen(word);  }

        word = strtok( NULL, " ,.;-()\n\r" );
    }

    // different criteria may be used for deciding whether the tried
    // key was a valid one. here we identify it as valid if words in
    // the decrypted message that can be located in the dictionary account
    // for more than half of the original message length.
    if (nmatches > len * 0.50)
        return 1;

    return 0;
}

void decrypt32( unsigned char *inp, uint32_t key, unsigned char *decoded ) {

    int i, iend, oend;
    uint32_t block;
    uint32_t a, b, c, d, magnitude, polarity, xor;

    srand(key);

    iend = 0;
    decoded[0] = 0; // C strings are zero-terminated
    oend = 0;

    /* main loop for decoding -- all 4 bytes are valid */
    while (     (a = inp[iend++]) != 0
            &&  (b = inp[iend++]) != 0
            &&  (c = inp[iend++]) != 0
            &&  (d = inp[iend++]) != 0
    ) {
        // printf("a = %x, b = %x, c = %x, d=%x\n", a, b, c, d);

        polarity = rand() % 2;
        magnitude = rand() % 32;

        block = (d << 24) | (c << 16) | (b << 8) | a;

        if (polarity)
            block = ROTR32(block, magnitude);
        else
            block = ROTL32(block, magnitude);

        xor =   (rand() % 256 << 24) | (rand() % 256 << 16)
              | (rand() % 256 << 8)  | (rand() % 256);

        block ^= xor;

        decoded[oend++] = block;
        decoded[oend++] = (block = block >> 8);
        decoded[oend++] = (block = block >> 8);
        decoded[oend++] = (block = block >> 8);
        decoded[oend] = 0;

        // printf("p = %d, mag = %d, xor = %d\n", polarity, magnitude, xor);
    }

    /* end cases */
    if ( a != 0 && b != 0 && c != 0 && d == 0 ) {

        polarity = rand() % 2;
        magnitude = rand() % 24;
        block = (c << 16) | (b << 8) | a;

        if (polarity)
            block = ROTR24(block, magnitude);
        else
            block = ROTL24(block, magnitude);

        xor = (rand() % 256 << 16) | (rand() % 256 << 8) | (rand() % 256);

        block ^= xor;

        decoded[oend++] = block;
        decoded[oend++] = (block = block >> 8);
        decoded[oend++] = (block = block >> 8);
        decoded[oend] = 0;

    } else if ( a != 0 && b != 0 && c == 0 ) {

        polarity = rand() % 2;
        magnitude = rand() % 16;
        block = (b << 8) | a;

        if (polarity)
            block = ROTR16(block, magnitude);
        else
            block = ROTL16(block, magnitude);

        xor = (rand() % 256 << 8) | (rand() % 256);

        block ^= xor;

        decoded[oend++] = block;
        decoded[oend++] = (block = block >> 8);
        decoded[oend] = 0;

    } else if ( a != 0 && b == 0 ) {

        polarity = rand() % 2;
        magnitude = rand() % 8;
        block = a;

        if (polarity)
            block = ROTR8(block, magnitude);
        else
            block = ROTL8(block, magnitude);

        xor = rand() % 256;
        block ^= xor;

        decoded[oend++] = block;
        decoded[oend] = 0;
    }
}


int main( int argc,char *argv[] ) {

    MPI_Init(&argc, &argv);
    
    char a;
    char outfilename[100];
    unsigned char encrypted[1000], decrypted[1000], dcopy[1000];
    FILE *fin, *fout;
    int success = 0; 
    uint32_t i, len; 
    double tstart, tend;
    int my_rank, procs, file_error = 0;  // file_error stores whether the input file is read successfully 
    
    MPI_Comm_size(MPI_COMM_WORLD, &procs);  // get the total number of processes 
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);  // get the rank of my process
    
    if (my_rank == 0)  // 0 process is responsible for reading file
    {
        printf("\n\nx0r 32-bit code breaker\n\n");
        
        if (argc == 1) {
            fprintf(stderr, "ERROR: No file(s) supplied.\n");
            fprintf(stderr, "USAGE: This program requires a filename to be provided as argument!");
            file_error = 1;
        }
        else
        {
            printf("decrypting file %s by trying all possible keys...\n", argv[1]);
            printf("To quit, press ctrl + c\n\n");
            printf("Status:\n");

            if ( (fin = fopen(argv[1], "r")) == NULL ) {
                fprintf(stderr, "ERROR: Could not open: %s\n", argv[1]);
                file_error = 1;
            }
            else
            {
                len = 0;
                while ( (a = fgetc(fin)) != EOF )
                    encrypted[len++] = a;
                encrypted[len] = 0;
                fclose(fin);
            }
        }
    }
    MPI_Bcast(&file_error, 1, MPI_INT, 0, MPI_COMM_WORLD);  // broadcast whether the input file is read successfully 
    if (file_error)  // file was not read, terminate all the processes 
    {
        MPI_Finalize();
        exit(1); 
    }
    
    // printf("encrypted: ");
    // for (i=0; i<end; ++i)
    //  printf("[%d]: %c", i, encrypted[i]);

    MPI_Barrier(MPI_COMM_WORLD);  // Synchoronize for time measurment 
    tstart = MPI_Wtime();
    
    MPI_Bcast(&len, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);  // broadcast the length of the string 
    MPI_Bcast(encrypted, len+1, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);  // broadcast the encrypted string 
    
    // calculate the range of trial keys 
    uint32_t my_start, my_end, task, limit; 
    limit = (uint32_t)(pow( 2, sizeof(uint32_t)*8 )); 
    task = limit / procs; 
    my_start = my_rank * task; 
    if (my_rank + 1 < procs) 
    {
        my_end = (my_rank + 1) * task;
    }
    else
    {
        my_end = limit; 
    }
    
    MPI_Request req;
    int flag = 0, dummy = 0; 
    MPI_Irecv(&dummy, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &req);  // receive whether other processes find the correct key
    
    // try to decrypt
    for ( i = my_start; i < my_end; ++i ) {

        decrypt32( encrypted, i, decrypted );

        // printf( "i = %u - decrypted: %s\n", i, decrypted );

        strcpy(dcopy, decrypted);

        if (isValid(dcopy, len)) {
            success = flag = 1; 
            
            for (int j = 0; j < procs; j++) 
               if (j != my_rank) MPI_Send(&success, 1, MPI_INT, j, 0, MPI_COMM_WORLD);  // tell my success to other processes
            
            break; 
        }
        
        MPI_Test(&req, &flag, MPI_STATUS_IGNORE);  // test whether other processes find the correct key 
        if (flag) break; 
    }
    
    MPI_Barrier(MPI_COMM_WORLD); 
    if (!flag) MPI_Test(&req, &flag, MPI_STATUS_IGNORE);  // final test of whether other processes find the key 
    
    tend = MPI_Wtime(); 

    if (success) { 
        // I find the correct key, terminate the program 
        
        printf("\nTime elapsed: %.6e seconds\n", tend - tstart); 
        
        sprintf(outfilename, "%s.out", argv[1]);
        fout = fopen(outfilename, "w");
        fprintf(fout, "%s", decrypted);
        printf("\nFile decrypted successfully using key %u\n", (unsigned) i);
        printf("See the file %s\n\n\n", outfilename);
        fclose(fout);
    
        return 0;
    }

    if (flag)
    {
        // other process finds the key, terminate the program  
        MPI_Finalize();
        return 0;
    }
    
    // cannot decrypt, terminate the program 
    if (my_rank == 0) 
    {
        printf("\nTime elapsed: %.6e seconds\n", tend - tstart); 
        printf("\nWARNING: File could not be decrypted.\n\n\n");
    }
    MPI_Finalize(); 
    return 1;
}
