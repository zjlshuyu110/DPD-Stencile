#include "quicksort.h"


int main(int argc, char const *argv[])
{
	return 0;
}


int check_and_print(NUMBER *elements, int n, char *file_name)
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


int read_input(char *file_name, NUMBER **elements)
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
