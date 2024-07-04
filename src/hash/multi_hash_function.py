from .hash_function import HashFunction

class MultiHashFunction:
    def __init__(self):
        self.hash_functions = []

    def add_hash_function(self, hash_function):
        self.hash_functions.append(hash_function)

    def hash(self, s):
        return [hf.hash(s) for hf in self.hash_functions]
