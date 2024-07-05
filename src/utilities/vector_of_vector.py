# vector_of_vector.py

def get_vector_of_vector(data, delimiter=','):
    result = []
    for line in data.splitlines():
        if line:
            result.append([int(x) for x in line.split(delimiter)])
    return result

def print_vector_of_vector(vector_of_vector):
    for vector in vector_of_vector:
        print(" ".join(map(str, vector)))
