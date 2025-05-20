import random

# Generate 100 random integers between 0 and 999
numbers = [str(random.randint(0, 999)) for _ in range(100)]

with open("input.txt", "w") as f:
    f.write(' '.join(numbers))