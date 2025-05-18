#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>

// ---------------- HELPER DEFINITIONS ----------------
#define MAX_SIZE 1000

// Swap for bubble sort
void swap(double *p, double *q) {
    double t = *p;
    *p = *q;
    *q = t;
}

// Bubble sort (demo only; replace with qsort or a real quicksort if you want)
void local_sort(double a[], int n) {
    for (int i = 0; i < n - 1; i++) {
        for (int j = 0; j < n - i - 1; j++) {
            if (a[j] > a[j + 1]) {
                swap(&a[j], &a[j + 1]);
            }
        }
    }
}

// Generate random array for testing
void generate_random_array(double arr[], int size, int range) {
    for (int i = 0; i < size; i++) {
        arr[i] = rand() % range;
    }
}

// Print array
void print_array(const double arr[], int size) {
    for (int i = 0; i < size; i++) {
        printf("%.0f ", arr[i]);
    }
    printf("\n");
}

// Partition data into chunk pointers, each chunk has size = total_size / num_chunks
// void partition_data(double* array, int total_size, int num_chunks,
//                     double** chunk_ptrs, int* chunk_sizes)
// {
//     int chunk_size = total_size / num_chunks;
//     for(int i = 0; i < num_chunks; i++){
//         chunk_ptrs[i] = &array[i * chunk_size];
//         chunk_sizes[i] = chunk_size;
//     }
// }


void partition_data(double* array, int total_size, int num_chunks,
    double** chunk_ptrs, int* chunk_sizes)
    {
    int chunk_size = total_size / num_chunks;
    for(int i = 0; i < num_chunks; i++){
    chunk_sizes[i] = chunk_size;
    // Actually allocate a new block for this chunk
    chunk_ptrs[i] = (double*) malloc(chunk_size * sizeof(double));
    // Copy the data from the big array into this chunk
    for(int j = 0; j < chunk_size; j++){
    chunk_ptrs[i][j] = array[i*chunk_size + j];
    }
    }
    // If there's a remainder, you'll want to handle that carefully, too
    }


// Find median by sorting local data & picking middle
double find_median(double *array, int size) {
    local_sort(array, size);
    if (size % 2 == 0) {
        return (array[size/2 - 1] + array[size/2]) / 2.0;
    } else {
        return array[size/2];
    }
}

// Sort an array of doubles in ascending order (for array of medians)
void sort_medians(double *arr, int n) {
    // re-use local_sort or do your own
    local_sort(arr, n);
}

// Select pivot from array of medians (Strategy 3 from lectures)
double pivot_select(double *meds, int n) {
    // Sort medians
    sort_medians(meds, n);
    // Take mean of middle two if even, or middle if odd
    if (n % 2 == 0)
        return (meds[n/2 - 1] + meds[n/2]) * 0.5;
    else
        return meds[n/2];
}

// ------------------ STEP 3.2: findsplit (partition local data by pivot) ------------------
// We build new arrays for smaller[] and bigger[]. Then the caller "merges" them.
void findsplit(double *local_data, int local_size, double pivot,
               double **smaller_part, int *smaller_count,
               double **bigger_part,  int *bigger_count)
{
    double *small_arr = (double*) malloc(local_size * sizeof(double));
    double *big_arr   = (double*) malloc(local_size * sizeof(double));

    if(!small_arr || !big_arr){
        printf("Malloc error in findsplit!\n");
        exit(1);
    }

    int sCount = 0, bCount = 0;
    for(int i = 0; i < local_size; i++){
        if(local_data[i] <= pivot) {
            small_arr[sCount++] = local_data[i];
        } else {
            big_arr[bCount++]   = local_data[i];
        }
    }

    *smaller_part = small_arr;
    *bigger_part  = big_arr;
    *smaller_count = sCount;
    *bigger_count  = bCount;
}

// ------------------ STEP 3.3: Merge data with partner ------------------
// We'll do two functions: MergeLower and MergeUpper
// "MergeLower" => chunk i keeps small_i+small_j, chunk j keeps big_i+big_j
void MergeLower(double **data_i, int *size_i,
                double *small_i, int count_small_i,
                double *big_i,   int count_big_i,
                double **data_j, int *size_j,
                double *small_j, int count_small_j,
                double *big_j,   int count_big_j)
{
    // free the old arrays
    free(*data_i);
    free(*data_j);

    // i: small_i + small_j
    int newSizeI = count_small_i + count_small_j;
    double *newDataI = (double*) malloc(newSizeI * sizeof(double));
    for(int k=0; k < count_small_i; k++){
        newDataI[k] = small_i[k];
    }
    for(int k=0; k < count_small_j; k++){
        newDataI[count_small_i + k] = small_j[k];
    }
    local_sort(newDataI, newSizeI); // final local sort
    *data_i = newDataI;
    *size_i = newSizeI;

    // j: big_i + big_j
    int newSizeJ = count_big_i + count_big_j;
    double *newDataJ = (double*) malloc(newSizeJ * sizeof(double));
    for(int k=0; k < count_big_i; k++){
        newDataJ[k] = big_i[k];
    }
    for(int k=0; k < count_big_j; k++){
        newDataJ[count_big_i + k] = big_j[k];
    }
    local_sort(newDataJ, newSizeJ);
    *data_j = newDataJ;
    *size_j = newSizeJ;
}

// "MergeUpper" => chunk i keeps big_i+big_j, chunk j keeps small_i+small_j
void MergeUpper(double **data_i, int *size_i,
                double *small_i, int count_small_i,
                double *big_i,   int count_big_i,
                double **data_j, int *size_j,
                double *small_j, int count_small_j,
                double *big_j,   int count_big_j)
{
    free(*data_i);
    free(*data_j);

    // i: big_i + big_j
    int newSizeI = count_big_i + count_big_j;
    double *newDataI = (double*) malloc(newSizeI * sizeof(double));
    for(int k=0; k < count_big_i; k++){
        newDataI[k] = big_i[k];
    }
    for(int k=0; k < count_big_j; k++){
        newDataI[count_big_i + k] = big_j[k];
    }
    local_sort(newDataI, newSizeI);
    *data_i = newDataI;
    *size_i = newSizeI;

    // j: small_i + small_j
    int newSizeJ = count_small_i + count_small_j;
    double *newDataJ = (double*) malloc(newSizeJ * sizeof(double));
    for(int k=0; k < count_small_i; k++){
        newDataJ[k] = small_i[k];
    }
    for(int k=0; k < count_small_j; k++){
        newDataJ[count_small_i + k] = small_j[k];
    }
    local_sort(newDataJ, newSizeJ);
    *data_j = newDataJ;
    *size_j = newSizeJ;
}

// ------------------ Recursive Parallel Global Sort ------------------
// We'll do a single function that uses OpenMP tasks for recursion.
// void parallel_global_sort(double **chunks, int *sizes, int group_size)
// {
//     if (group_size <= 1) return;

//     // 3.1: pivot selection
//     //  - each chunk finds median in parallel
//     double *medians = (double*) malloc(group_size * sizeof(double));

//     // Create a task for each chunk to find its median
//     #pragma omp for
//     for(int i = 0; i < group_size; i++){
//         medians[i] = find_median(chunks[i], sizes[i]);
//     }

//     // Let a single thread select the pivot
//     double pivot;
// #pragma omp single
//     {
//         pivot = pivot_select(medians, group_size);
//     }

//     // barrier to ensure pivot is known by all
// #pragma omp barrier

//     // 3.2: findsplit in parallel
//     // We'll store smaller/bigger arrays in parallel arrays for each chunk
//     double *smallers[128], *biggers[128];  // assume group_size <= 128 for demo
//     int count_smallers[128], count_biggers[128];

// #pragma omp for
//     for(int i = 0; i < group_size; i++){
//         findsplit(chunks[i], sizes[i], pivot,
//                   &smallers[i], &count_smallers[i],
//                   &biggers[i], &count_biggers[i]);
//     }

//     // barrier so that all splits are done before merges
// #pragma omp barrier

//     // 3.3: merge in pairs: i <-> i+half
//     int half = group_size / 2;

//     // We'll do merges in parallel, so each pair i, i+half merges
// #pragma omp for
//     for(int i = 0; i < half; i++){
//         int partner = i + half;
//         // We keep the "lower" portion for chunk i => MergeLower
//         MergeLower(&chunks[i], &sizes[i],
//                    smallers[i], count_smallers[i],
//                    biggers[i], count_biggers[i],
//                    &chunks[partner], &sizes[partner],
//                    smallers[partner], count_smallers[partner],
//                    biggers[partner], count_biggers[partner]);
//     }

//     // barrier to ensure merges are done before recursion
// #pragma omp barrier

//     // free the smallers/biggers arrays
// #pragma omp for
//     for(int i = 0; i < group_size; i++){
//         free(smallers[i]);
//         free(biggers[i]);
//     }

//     free(medians);

//     // 3.4: recursion on left half and right half
//     // We'll spawn tasks to do these subgroups in parallel
//     #pragma omp single
//     {
//         printf("Creating tasks for left and right subgroups, group_size=%d\n", group_size);
//         #pragma omp task
//         {
//             printf("Starting left recursion group_size=%d\n", half);
//             parallel_global_sort(chunks, sizes, half);
//             printf("Finished left recursion group_size=%d\n", half);
//         }
//         #pragma omp task
//         {
//             printf("Starting right recursion group_size=%d\n", group_size - half);
//             parallel_global_sort(&chunks[half], &sizes[half], group_size - half);
//             printf("Finished right recursion group_size=%d\n", group_size - half);
//         }
//         #pragma omp taskwait
//         printf("All tasks done for group_size=%d\n", group_size);
//     }
//     // end single

//     // barrier to ensure recursion tasks are done
// #pragma omp barrier
// }


void parallel_global_sort(double **chunks, int *sizes, int group_size)
{
    // Base case
    if (group_size <= 1) return;

    // 3.1: pivot selection
    //  - each chunk finds its median in parallel tasks
    double *medians = (double*) malloc(group_size * sizeof(double));

    // Create a task for each chunk to find its median
    #pragma omp taskgroup
    {
        for (int i = 0; i < group_size; i++){
            #pragma omp task firstprivate(i) shared(chunks, sizes, medians)
            {
                medians[i] = find_median(chunks[i], sizes[i]);
            }
        }
    }

    // Let a single thread select the pivot from all medians
    double pivot;
    #pragma omp single
    {
        pivot = pivot_select(medians, group_size);
    }

    // 3.2: findsplit in parallel for each chunk
    double *smallers[128], *biggers[128]; 
    int count_smallers[128], count_biggers[128];

    #pragma omp taskgroup
    {
        for(int i = 0; i < group_size; i++){
            #pragma omp task firstprivate(i) shared(chunks, sizes, pivot, smallers, biggers, count_smallers, count_biggers)
            {
                findsplit(chunks[i], sizes[i], pivot,
                          &smallers[i], &count_smallers[i],
                          &biggers[i], &count_biggers[i]);
            }
        }
    }

    // 3.3: pairwise merge (i, i+half) => MergeLower
    int half = group_size / 2; // integer division
    int leftover = group_size - half;

    // Merge pairs in parallel tasks
    #pragma omp taskgroup
    {
        for(int i = 0; i < half; i++){
            int partner = i + half; 
            #pragma omp task firstprivate(i, partner) \
                            shared(chunks, sizes, smallers, biggers, count_smallers, count_biggers)
            {
                MergeLower(&chunks[i], &sizes[i],
                           smallers[i], count_smallers[i],
                           biggers[i], count_biggers[i],
                           &chunks[partner], &sizes[partner],
                           smallers[partner], count_smallers[partner],
                           biggers[partner], count_biggers[partner]);
            }
        }
    }

    // free all smallers/biggers arrays
    for(int i = 0; i < group_size; i++){
        free(smallers[i]);
        free(biggers[i]);
    }
    free(medians);

    // 3.4: recursion on left half and right half, each in its own task
    #pragma omp taskgroup
    {
        #pragma omp task
        {
            parallel_global_sort(chunks, sizes, half);
        }
        #pragma omp task
        {
            parallel_global_sort(&chunks[half], &sizes[half], leftover);
        }
    }
}


// ------------------ MAIN FUNCTION ------------------
// int main(int argc, char* argv[])
// {
//     // Command-line argument for threads
//     int num_threads = 8; // default
//     if (argc > 1) {
//         num_threads = atoi(argv[1]);
//     }
//     omp_set_num_threads(num_threads);

//     srand((unsigned)time(NULL));

//     // Generate data
//     int total_size = MAX_SIZE;
//     double *array = (double*) malloc(sizeof(double) * total_size);
//     generate_random_array(array, total_size, 1000);

//     printf("Initial array:\n");
//     print_array(array, total_size);

//     // Partition into chunks for "NUM_THREADS" or something.
//     // But let's keep it at some chunk_count = num_threads for a test.
//     int chunk_count = num_threads;
//     double* chunk_ptrs[256];  // up to 256 for demo
//     int chunk_sizes[256];

//     partition_data(array, total_size, chunk_count, chunk_ptrs, chunk_sizes);

//     for(int i=0; i<chunk_count; i++){
//         local_sort(chunk_ptrs[i], chunk_sizes[i]);
//     }

//     // Parallel region: we run parallel_global_sort in a single parallel session.
//     double start = omp_get_wtime();

// #pragma omp parallel
//     {
// #pragma omp single
//         {
//             parallel_global_sort(chunk_ptrs, chunk_sizes, chunk_count);
//         }
//     }

//     double end = omp_get_wtime();

//     // Print final chunk data
//     printf("\nChunks after parallel_global_sort:\n");
//     for(int i=0; i<chunk_count; i++){
//         print_array(chunk_ptrs[i], chunk_sizes[i]);
//     }

//     printf("\nTotal time: %f seconds using %d threads\n", (end - start), num_threads);

//     free(array);
//     return 0;
// }


int main(int argc, char* argv[])
{
    int num_threads = 4; // default
    if (argc > 1) {
        num_threads = atoi(argv[1]);
    }
    omp_set_num_threads(num_threads);

    srand((unsigned)time(NULL));

    // Generate data
    int total_size = MAX_SIZE;
    double *array = (double*) malloc(sizeof(double) * total_size);
    generate_random_array(array, total_size, 1000);

    printf("Initial array:\n");
    print_array(array, total_size);

    // Partition into chunks (1 chunk per thread for a simple test)
    int chunk_count = num_threads;
    double* chunk_ptrs[256];
    int chunk_sizes[256];

    partition_data(array, total_size, chunk_count, chunk_ptrs, chunk_sizes);

    //  local sort each chunk
    for(int i = 0; i < chunk_count; i++){
        local_sort(chunk_ptrs[i], chunk_sizes[i]);
    }

    double start = omp_get_wtime();

    #pragma omp parallel
    {
        // We enter parallel, but only one thread calls parallel_global_sort
        #pragma omp single
        {
            parallel_global_sort(chunk_ptrs, chunk_sizes, chunk_count);
        }
    }

    double end = omp_get_wtime();

    // Print final chunk data
    printf("\nChunks after parallel_global_sort:\n");
    for(int i=0; i<chunk_count; i++){
        print_array(chunk_ptrs[i], chunk_sizes[i]);
    }
    printf("\nTotal time: %f seconds using %d threads\n", (end - start), num_threads);

    free(array);
    return 0;
}