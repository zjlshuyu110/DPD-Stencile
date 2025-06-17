import random

total = 1_000_000
# Generate 100 random integers between 0 and 999
numbers = [str(random.randint(0, total)) for _ in range(total)]

with open("A3/input.txt", "w") as f:
    f.write(f"{total}" + " ")
    f.write(' '.join(numbers))