# utilities.py

import os

class Utilities:
    @staticmethod
    def parse_line(line, delimiter='/'):
        return line.strip().split(delimiter)

    @staticmethod
    def create_directory(path):
        if not os.path.exists(path):
            os.makedirs(path)

    @staticmethod
    def file_exists(filepath):
        return os.path.isfile(filepath)

    @staticmethod
    def read_file(filepath):
        with open(filepath, 'r') as file:
            return file.read()

    @staticmethod
    def write_file(filepath, data):
        with open(filepath, 'w') as file:
            file.write(data)

def print_vector(vector):
    print(" ".join(map(str, vector)))

def print_vector_of_vector(vector_of_vector):
    for vector in vector_of_vector:
        print_vector(vector)
