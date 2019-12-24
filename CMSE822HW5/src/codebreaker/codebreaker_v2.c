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
#define ROTL8(x,y) ((x<<y)|(x>>(8-y)))
#define ROTR8(x,y) ((x>>y)|(x<<(8-y)))

int num_entries = 8;
const char *dictionary[] = {"bird", "campus", "class", "of", 
    "spring", "sun", "the", "tree"}; 


// checks if a significant portion of the words resulting from 
// the decryption with the tried key are found in the dictionary 
int isValid(char *decoded, int len)
{
    int i, flag;
    int nmatches = 0;
    char *word;

    word = strtok(decoded," ,.;-()\n\r");
    while (word != NULL) { 
        flag = 0;
        for (i=0; i < num_entries; ++i) {
            if(strcmp(word, dictionary[i]) == 0) {
                flag = 1;
                break;
            }	
        }
        if (flag)
        {
            nmatches += strlen(word);
            //printf("MATCH: %s (%d)\n",word,strlen(word));
        }
        word = strtok(NULL," ,.;-()\n\r");
    }

    // different criteria may be used for deciding whether the tried 
    // key was a valid one. here we identify it as valid if words in 
    // the decrypted message that can be located in the dictionary account 
    // for more than half of the original message length.
    if (nmatches > len * 0.50)
        return 1;

    return 0;
}


void decrypt32(unsigned char *inp, unsigned key, unsigned char *decoded, unsigned len)
{
    int iend, oend;
    unsigned block;
    unsigned a, b, c, d, magnitude, polarity, xor;
    unsigned remaining;

    srand(key);

    remaining = len;
    iend = 0;
    oend = 0;        

    /* main loop for decoding -- all 4 bytes are valid */
    while(remaining >= 4) {
        a = inp[iend++];
        b = inp[iend++];
        c = inp[iend++];
        d = inp[iend++];

        //printf("a = %x, b = %x, c = %x, d=%x\n", a, b, c, d);
        polarity = rand()%2;
        magnitude = rand()%32;
        block = ((d<<24) | (c<<16) |(b<<8) | a);

        if (polarity) 
            block = ROTR32(block,magnitude);
        else block = ROTL32(block,magnitude);

        xor = ((rand()%256<<24) | (rand()%256<<16) | 
                (rand()%256<<8) | rand()%256);
        block ^= xor;

        decoded[oend++] = block;
        decoded[oend++] = (block=block>>8);
        decoded[oend++] = (block=block>>8);
        decoded[oend++] = (block=block>>8);
        //printf("p = %d, mag = %d, xor = %d\n", polarity, magnitude, xor);

        remaining -= 4;
    }

    /* end cases */
    if (remaining == 3) {
        a = inp[iend++];
        b = inp[iend++];
        c = inp[iend++];

        polarity = rand()%2;
        magnitude = rand()%24;
        block = ((c<<16) |(b<<8) | a);

        if (polarity) 
            block = ROTR24(block,magnitude);
        else block = ROTL24(block,magnitude);

        xor = ((rand()%256<<16) | (rand()%256<<8) | rand()%256);
        block ^= xor;

        decoded[oend++] = block;
        decoded[oend++] = (block=block>>8);
        decoded[oend++] = (block=block>>8);
    }
    else if (remaining == 2) {
        a = inp[iend++];
        b = inp[iend++];

        polarity = rand()%2;
        magnitude = rand()%16;
        block = ((b<<8) | a);

        if (polarity) 
            block = ROTR16(block,magnitude);
        else block = ROTL16(block,magnitude);

        xor = ((rand()%256<<8) | rand()%256);
        block ^= xor;

        decoded[oend++] = block;
        decoded[oend++] = (block=block>>8);
    }
    else if (remaining == 1) {
        a = inp[iend++];

        polarity = rand()%2;
        magnitude = rand()%8;
        block = (a);

        if (polarity) 
            block = ROTR8(block,magnitude);
        else block = ROTL8(block,magnitude);

        xor = (rand()%256);
        block ^= xor;

        decoded[oend++] = block;
    }
}


void decrypt8(unsigned char *inp, unsigned key, unsigned char *decoded, unsigned len)
{
    int iend, oend;
    unsigned char a, magnitude, block, xor;
    unsigned remaining;

    remaining = len;
    iend = 0;
    decoded[0] = 0; // C strings are zero-terminated
    oend = 0;        

    srand(key);

    /* main loop for encoding */
    while(remaining > 0) {
        a = inp[iend++];
        magnitude = rand()%8;
        xor = rand()%256;
        block = a;
        block = ROTR8(block, magnitude);
        block ^= xor;
        decoded[oend++] = block;
        decoded[oend] = 0;
        //    printf("%x --> %x, mag: %d, xor: %x\n", a, block, magnitude, xor);

        --remaining;
    }

}


int main(int argc,char *argv[])
{
    MPI_Init(&argc, &argv);
    
    char a;
    char outfilename[100];
    unsigned char encrypted[1000], decrypted[1000], dcopy[1000];
    FILE *fin, *fout;
    int success = 0;
    unsigned i, len;
    double tstart, tend;
    int my_rank, procs;
    
    MPI_Comm_size(MPI_COMM_WORLD, &procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    
    if (my_rank == 0)
    {
        printf("\n\nx0r 32-bit code breaker\n\n");

        if (argc == 1) {
            fprintf(stderr, "ERROR: No file(s) supplied.\n");
            fprintf(stderr, "USAGE: This program requires a filename \
                    to be provided as argument!");
            exit(1);
        }

        printf("decrypting file %s by trying all possible keys...\n", argv[1]);
        printf("To quit, press ctrl + c\n\n");
        printf("Status:\n");

        if ((fin = fopen(argv[1], "rb")) == NULL) {
            fprintf(stderr, "ERROR: Could not open: %s\n", argv[1]);
            exit(1);
        }
        
        fseek(fin, 0, SEEK_END);
        len = ftell(fin);
        fseek(fin, 0, SEEK_SET);
        fread(encrypted, len, 1, fin);
        printf("[LOG] %u encrypted bytes read.\n", len);
        fclose(fin);
    }

    //  printf("encrypted: "); 
    //  for (i=0; i<end; ++i)
    //    printf("[%d]: %c", i, encrypted[i]);
    
    tstart = MPI_Wtime();
    
    MPI_Bcast(&len, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(encrypted, len+1, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
    /*
    unsigned my_start, my_end, task, limit;
    limit = (unsigned)(pow(2,sizeof(unsigned)*8));
    task = limit/procs;
    my_start = my_rank*task;
    if (my_rank+1<procs)
    {
        my_end = (my_rank+1)*task;
    }
    else
    {
        my_end = limit;
    }
    */
    MPI_Request req;
    int flag = 0;
    int dummy = 0;
    MPI_Irecv(&dummy, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &req);
    
    for (i=my_rank; i<pow(2,sizeof(unsigned)*8); i=i+procs) {
        decrypt32(encrypted, i, decrypted, len);
        //    printf("i=%d - decrypted: %s\n", i, decrypted); 

        memcpy(dcopy, decrypted, len * sizeof(unsigned char));

        if (isValid(dcopy, len)) {
            success = 1;
            flag = 1;
            for (int j = 0; j<procs; j++)
                if (j != my_rank)
                    MPI_Send(&success, 1, MPI_INT, j, 0, MPI_COMM_WORLD);
            break;
        }
        MPI_Test(&req, &flag, MPI_STATUS_IGNORE);
        if (flag) break;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (!flag) MPI_Test(&req, &flag, MPI_STATUS_IGNORE);
    
    tend = MPI_Wtime();
    
    
    if (success) {
        printf("\nTime elapsed: %.6e seconds\n", tend - tstart);
        sprintf(outfilename, "%s.out", argv[1]);
        fout = fopen(outfilename, "w");
        fprintf(fout, "%s", decrypted); 
        printf("\nFile decrypted successfully using key %u\n", i);
        printf("See the file %s\n\n\n", outfilename);
        fclose(fout);
        MPI_Finalize();
        return 0;
    }
    
    if (flag)
    {
        // other process finds the key, terminate the program
        MPI_Finalize();
        return 0;
    }
    
    
    if (my_rank == 0)
    {
        printf("\nTime elapsed: %.6e seconds\n", tend - tstart);
        printf("\nWARNING: File could not be decrypted.\n\n\n");
    }
    MPI_Finalize();
    return 1;
}
