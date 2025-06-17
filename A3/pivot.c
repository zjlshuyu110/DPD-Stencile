#include "pivot.h"
#include <stdlib.h>

// Compares two integers for qsort
int compare(const void *v1, const void *v2) {
    int a = *(const int *)v1;
    int b = *(const int *)v2;
    return a - b;
}

// Returns the index of the first element larger than val, or n if none are found
int get_larger_index(int *elements, int n, int val) {
    for (int i = 0; i < n; ++i) {
        if (elements[i] > val)
            return i;
    }
    return n;
}

// Returns the median value from the sorted elements
int get_median(int *elements, int n) {
    if (n == 0) return 0; // or handle error
    if (n % 2 == 1)
        return elements[n / 2];
    else
        return (elements[n / 2 - 1] + elements[n / 2]) / 2;
}

// Selects a pivot based on the strategy and elements provided
int select_pivot(int pivot_strategy, int *elements, int n, MPI_Comm communicator) {
    switch (pivot_strategy) {
        case MEDIAN_ROOT:
            return select_pivot_median_root(elements, n, communicator);
        case MEAN_MEDIAN:
            return select_pivot_mean_median(elements, n, communicator);
        case MEDIAN_MEDIAN:
            return select_pivot_median_median(elements, n, communicator);
        default:
            return select_pivot_smallest_root(elements, n, communicator);
    }
}

// Selects the pivot using the median-of-medians strategy
int select_pivot_median_root(int *elements, int n, MPI_Comm communicator) {
	int rank;
	MPI_Comm_rank(communicator, &rank);
	int median = 0;
	if (rank == 0) median = get_median(elements, n);
	MPI_Bcast(&median, 1, MPI_INT, 0, communicator);

	return get_larger_index(elements, n, median);
}

// Selects the pivot using the mean-of-medians strategy (placeholder)
int select_pivot_mean_median(int *elements, int n, MPI_Comm communicator) {
	// Placeholder: implement mean of medians logic as needed
	int rank, size;
	MPI_Comm_rank(communicator, &rank);
	MPI_Comm_size(communicator, &size);
	int local_median = get_median(elements, n);
	int *medians;
	if (rank == 0) {
		medians = malloc(size * sizeof(int));
		if (medians == NULL) {
			printf("Error allocating buffers for medians\n");
			MPI_Abort(MPI_COMM_WORLD, 2);
		}
	}

	// Gather all local median to medians
	MPI_Gather(&local_median, 1, MPI_INT, medians, 1, MPI_INT, 0, communicator);
	int mean = 0;
	if (rank == 0) {
		int sum = 0;
		for (int i=0; i<size; i++)
			sum += medians[i];
		mean = sum / size;
		free(medians);
	}
	MPI_Bcast(&mean, 1, MPI_INT, 0, communicator);
	return get_larger_index(elements, n, mean);
}

static int compare_qsort(const void *a, const void *b)
{
	return (*(int*)a - *(int*)b);
}

// Selects the pivot using the median-of-medians strategy (placeholder)
int select_pivot_median_median(int *elements, int n, MPI_Comm communicator) {
	// Placeholder: implement median of medians logic as needed
	int rank, size;
	MPI_Comm_rank(communicator, &rank);
	MPI_Comm_size(communicator, &size);
	int local_median = get_median(elements, n);
	int *medians;
	if (rank == 0) {
		medians = malloc(size * sizeof(int));
		if (medians == NULL) {
			printf("Error allocating buffers for medians\n");
			MPI_Abort(MPI_COMM_WORLD, 2);
		}
	}

	// Gather all local median to medians
	MPI_Gather(&local_median, 1, MPI_INT, medians, 1, MPI_INT, 0, communicator);
	int median = 0;
	if (rank == 0) {
		qsort(medians, size, sizeof(int), compare_qsort);
		median = get_median(medians, size);
		free(medians);
	}
	MPI_Bcast(&median, 1, MPI_INT, 0, communicator);
	return get_larger_index(elements, n, median);
}

// Selects the pivot as the smallest element's index (fallback strategy)
int select_pivot_smallest_root(int *elements, int n, MPI_Comm communicator) {
    int rank;
	MPI_Comm_rank(communicator, &rank);
	int smallest = 0;
	if (rank == 0) smallest = elements[0];
	MPI_Bcast(&smallest, 1, MPI_INT, 0, communicator);
	return get_larger_index(elements, n, smallest);
}