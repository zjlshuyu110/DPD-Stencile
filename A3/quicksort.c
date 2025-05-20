#include "quicksort.h"


int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);

	if (argc != 4) {
		fprintf(stderr, "quicksort <input_file_name> <output_file_name> <pivot_strategy 1, 2, 3>");
		MPI_Abort(MPI_COMM_WORLD, 2);
	}
	const char* input_file_name = argv[1];
	if (input_file_name == NULL) {
		fprintf(stderr, "Please provide the input file name!\n");
		MPI_Abort(MPI_COMM_WORLD, 2);
	}
	const char* output_file_name = argv[2];
	if (output_file_name == NULL) {
		fprintf(stderr, "Please provide the output file name!\n");
		MPI_Abort(MPI_COMM_WORLD, 2);
	}
	const int pivot_strategy = atoi(argv[3]);
	if (pivot_strategy != 1 && pivot_strategy != 2 && pivot_strategy != 3) {
		fprintf(stderr, "Pivot Strategy must be choosen among 1, 2, or 3!\n");
		MPI_Abort(MPI_COMM_WORLD, 2);
	}

	int* global_elements;
	int* elements;
	int n = 0;

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank == 0) n = read_input(input_file_name, &global_elements);

	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

	int local_n = distribute_from_root(global_elements, n, &elements);

	serial_sort(elements, local_n);

	local_n = global_sort(&elements, local_n, MPI_COMM_WORLD, pivot_strategy);

	gather_on_root(global_elements, elements, local_n);

	if (rank == 0) {
		check_and_print(global_elements, n, output_file_name);
		free(global_elements);
	}

	free(elements);

	MPI_Finalize();
	return 0;
}


int check_and_print(int *elements, int n, const char *file_name)
{
	if (!sorted_ascending(elements, n)) {
		printf("Error! Array failed to sort!\n");
		return -1;
	}

	FILE *file = fopen(file_name, "w");  // "w" for text mode
	if (file == NULL) {
		printf("Error! Can't open the file %s\n", file_name);
		return -2;
	}

	// Write each element to the file
	for (int i = 0; i < n; ++i) {
		if (fprintf(file, "%d ", elements[i]) < 0) {
			fprintf(stderr, "Failed to write element %d to file.\n", i);
			fclose(file);
			return -2;
		}
	}

	fclose(file);
	return 0;
}

int global_sort(int **elements, int n, MPI_Comm comm, int pivot_strategy)
{
	int rank, size;
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &size);

	// Require even number of processes for pairing
	if (size > 1 && size % 2 != 0) {
		if (rank == 0) fprintf(stderr, "global_sort: Number of processes must be even. \n");
		MPI_Abort(comm, 1);
	}

	// Base case: single process
	if (size == 1) {
		return n;
	}

	// 3.1 Select pivot on rank 0 and broadcast
	int pivot;
	if (rank == 0) {
		int pivot_index = select_pivot(pivot_strategy, *elements, n, comm);
		pivot = (*elements)[pivot_index];
	}
	MPI_Bcast(&pivot, 1, MPI_INT, 0, comm);

	// 3.2 Partition locally around pivot.
	// Use B-search to find the split point of the array.
	int left = 0, right = n;
	int mid = 0;
	while (left < right) {
		mid = left + (right - left) / 2;
		if ((*elements)[mid] < pivot)
			left = mid + 1;
		else
			right = mid;
	}
	const int left_arr_size  = left;
	const int right_arr_size = n - left;

	// 3.3 Split processes into two groups
	int color = (rank < size/2) ? 0 : 1;
	MPI_Comm sub_comm;
	MPI_Comm_split(comm, color, rank, &sub_comm);

	int partner_rank = (color == 0) ? rank + size/2 : rank - size/2;

	// Exchange sizes.
	int send_count = (color == 0) ? right_arr_size : left_arr_size;
	int recv_count;
	MPI_Sendrecv(&send_count, 1, MPI_INT, partner_rank, 0,
					&recv_count, 1, MPI_INT, partner_rank, 0,
					comm, MPI_STATUS_IGNORE);

	// Exchange data.
	int *recv_buff = malloc(recv_count * sizeof(int));
	if (recv_buff == NULL) {
		printf("Malloc failed on rank %d\n", rank);
		MPI_Abort(MPI_COMM_WORLD, 2);
	}
	const int *send_ptr = (color == 0) ? (*elements + left_arr_size) : *elements;
	MPI_Sendrecv(send_ptr, send_count, MPI_INT, partner_rank, 0,
					recv_buff, recv_count, MPI_INT, partner_rank, 0,
					comm, MPI_STATUS_IGNORE);

	// 3.4 Merge the two sorted runs into one sorted array.
	int new_n = ((color == 0) ? left_arr_size : right_arr_size) + recv_count;
	int *new_elements = malloc(new_n * sizeof(int));
	if (new_elements == NULL) {
		printf("Allocation of new elements failed on rank %d\n", rank);
		MPI_Abort(MPI_COMM_WORLD, 2);
	}

	if (color == 0) {
		merge_ascending(*elements, left_arr_size,
						recv_buff, recv_count,
						new_elements);
	} else {
		merge_ascending(recv_buff, recv_count,
						*elements + left_arr_size, right_arr_size,
						new_elements);
	}

	free(*elements);
	free(recv_buff);
	*elements = new_elements;

	// 4 Recursive call.
	int final_n = global_sort(elements, new_n, sub_comm, pivot_strategy);
	MPI_Comm_free(&sub_comm);

	return final_n;
}

void merge_ascending(int *v1, int n1, int *v2, int n2, int *result)
{
	int i = 0, j = 0, k = 0;
	while (i < n1 && j < n2) {
		result[k++] = (v1[i] <= v2[j]) ? v1[i++] : v2[j++];
	}
	while (i < n1) result[k++] = v1[i++];
	while (j < n2) result[k++] = v2[j++];
}

int read_input(const char *file_name, int **elements)
{
	FILE *file = fopen(file_name, "r");
	if (file == NULL) {
		fprintf(stderr, "Unable to read the input file %s!\n", file_name);
		return -2;
	}

	int n;
	if (fscanf(file, "%d", &n) != 1 || n < 1) {
		fprintf(stderr, "Failed to read array size\n");
		fclose(file);
		return -2;
	}

	*elements = malloc(n * sizeof(NUMBER));
	if (!*elements) {
		fprintf(stderr, "Error allocating memory for %d elements!\n", n);
		fclose(file);
		return -2;
	}

	for (int i = 0; i < n; i++) {
		if (fscanf(file, "%d", &(*elements)[i]) != 1) {
			fprintf(stderr, "Failed to read element %d\n", i);
			free(*elements);
			fclose(file);
			return -2;
		}
	}

	fclose(file);

	return n;
}

int sorted_ascending(int *elements, int n)
{
	for (int i=1; i<n; i++) {
		if (elements[i] < elements[i-1])
			return 0;
	}
	return 1;
}

void swap(int *e1, int *e2)
{
	int temp = *e1;
	*e1 = *e2;
	*e2 = temp;
}

void serial_sort(int *elements, int n)
{
	for (int i = 1; i < n; i++) {
		int key = elements[i];
		int j = i - 1;

		while (j >= 0 && elements[j] > key) {
			elements[j + 1] = elements[j];
			j--;
		}
		elements[j + 1] = key;
	}
}

int distribute_from_root(int *all_elements, int n, int **my_elements)
{   // Distribute elements from root to all processes
    //to do: check if all_elements is NULL
    // all_elements is the buffer on root

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Calculate counts and displacements for each process
    int *counts = malloc(size * sizeof(int));
    int *displs = malloc(size * sizeof(int));
    int base = n / size;
    int rem = n % size;
    for (int i = 0, disp = 0; i < size; ++i) {
        counts[i] = base + (i < rem ? 1 : 0);
        displs[i] = disp;
        disp += counts[i];
    }

    int local_n = counts[rank];
    *my_elements = malloc(local_n * sizeof(int));

    MPI_Scatterv(
        all_elements, counts, displs, MPI_INT,
        *my_elements, local_n, MPI_INT,
        0, MPI_COMM_WORLD
    );

    free(counts);
    free(displs);
    return local_n;
}

void gather_on_root(int *all_elements, int *my_elements, int local_n)
{
    // Gather all elements on root
    // all_elements is the buffer on root
    // my_elements is the local array on each process
    // local_n is the number of elements in my_elements
    // all_elements will be allocated on root
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Gather counts from all processes
    int *counts = NULL;
    int *displs = NULL;
    if (rank == 0) {
        counts = malloc(size * sizeof(int));
    }
    MPI_Gather(&local_n, 1, MPI_INT, counts, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        displs = malloc(size * sizeof(int));
        displs[0] = 0;
        for (int i = 1; i < size; ++i) {
            displs[i] = displs[i-1] + counts[i-1];
        }
    }

    MPI_Gatherv(
        my_elements, local_n, MPI_INT,
        all_elements, counts, displs, MPI_INT,
        0, MPI_COMM_WORLD
    );

    if (rank == 0) {
        free(counts);
        free(displs);
    }
}
