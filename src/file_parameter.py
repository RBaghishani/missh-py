class FileParameter:
    def __init__(self):
        self.spaced_qmers = []

    def init(self, file1, file2=""):
        # Initialize the parameters based on the files
        # For now, return True as in C++ code
        return True

    def add_spaced_qmer(self, qmer1, qmer2):
        self.spaced_qmers.append((qmer1, qmer2))
