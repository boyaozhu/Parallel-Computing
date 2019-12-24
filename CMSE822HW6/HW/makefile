all: bucket_sort_v1 bucket_sort_v2 bucket_sort_v3

bucket_sort_v1: bucket_sort_v1.c
	mpicc bucket_sort_v1.c -o bucket_sort_v1 -O3

bucket_sort_v2: bucket_sort_v2.c
	mpicc bucket_sort_v2.c -o bucket_sort_v2 -O3

bucket_sort_v3: bucket_sort_v3.c
	mpicc -lm bucket_sort_v3.c -o bucket_sort_v3 -O3


clean: 
	rm bucket_sort_v1 bucket_sort_v2 bucket_sort_v3 slurm* shab* core*