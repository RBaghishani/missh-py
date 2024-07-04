from enum import Enum

class HashType(Enum):
    MD5 = 1
    SHA1 = 2
    SHA256 = 3

def hash_type_to_string(hash_type):
    if hash_type == HashType.MD5:
        return "MD5"
    elif hash_type == HashType.SHA1:
        return "SHA1"
    elif hash_type == HashType.SHA256:
        return "SHA256"
    else:
        return "Unknown"
