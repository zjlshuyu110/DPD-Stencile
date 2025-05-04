#include "quicksort.h"


int main(int argc, char const *argv[])
{
	const char* input_file_name = "test.txt";
	const char* output_file_name = "result.txt";
	int* global_elements;
	int n = 0;

	MPI_Init(&argc, &argv);

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank == 0) n = read_input(input_file_name, &global_elements);
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

	global_sort(&global_elements, n, MPI_COMM_WORLD, 1);

	if (rank == 0) check_and_print(global_elements, n, output_file_name);

	return 0;
}


int check_and_print(int *elements, int n, char *file_name)
{
	if (!sorted_ascending(elements, n)) {
		printf("Error! Array failed to sort!\n");
		return -1;
	}

	FILE *file = fopen(file_name, "wb");
	if (file == NULL) {
		printf("Error! Can't open the file %s\n", file_name);
		return -2;
	}

	if (fwrite(&n, sizeof(int), 1, file) != 1) {
		fprintf(stderr, "Failed to write array size to file.\n");
		fclose(file);
		return -2;
	}

	if (fwrite(elements, sizeof(NUMBER), n, file) != (size_t)n) {
		fprintf(stderr, "Failed to write array data to file.\n");
		fclose(file);
		return -2;
	}

	fclose(file);
}


int global_sort(int **elements, int n, MPI_Comm, int pivot_strategy)
{
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	int *local_elements;
	int local_n = distribute_from_root(*elements, n, &local_elements);
	if (local_n < 1) {
		return -2;
	}

	int pivot_index;
	if (rank == 0) pivot_index = select_pivot(pivot_strategy, *local_elements, local_n, MPI_COMM_WORLD);
	int pivot = local_elements[pivot_index];
	MPI_Bcast(&pivot, 1, MPI_INT, 0, MPI_COMM_WORLD);



	gather_on_root(elements, local_elements, local_n);
	free(local_elements);
	return 0;
}


int read_input(char *file_name, int **elements)
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

void swap(int *e1, int *e2)
{
	int temp = *e1;
	*e1 = *e2;
	*e2 = temp;
}
