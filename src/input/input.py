import os

class FileType:
    Fasta, Fastq, UnknownFile = range(3)

class PairType:
    SingleEnd, PairedEnd, UnknownPair = range(3)

class SingleEndFile:
    def __init__(self):
        self.reset()

    def reset(self):
        self.path = ""
        self.path_parse = []
        self.file_type = FileType.UnknownFile
        self.read_delimiter = '>'
        self.open_correct = False

    def init(self, path):
        self.reset()

        if not os.path.isfile(path):
            print(f"Fail to open: {path}")
            return False

        with open(path) as file:
            line = file.readline().strip()
            if not line:
                return False

            self.path = path
            self.path_parse = path.split("/")
            self.open_correct = True

            if line[0] == '>':
                self.file_type = FileType.Fasta
                self.read_delimiter = '>'
            elif line[0] == '@':
                self.file_type = FileType.Fastq
                self.read_delimiter = '@'
            else:
                return False
            return True

    def is_correct(self):
        return self.open_correct and self.file_type != FileType.UnknownFile

    def get_path(self):
        return self.path

    def get_path_parse(self):
        return self.path_parse

    def get_directory(self):
        return "/".join(self.path_parse[:-1]) + "/"

    def get_filename(self):
        return self.path_parse[-1]

    def get_filename_without_ext(self):
        filename_without_ext = self.path_parse[-1]
        pos = filename_without_ext.rfind(".")
        if pos != -1:
            filename_without_ext = filename_without_ext[:pos]
        return filename_without_ext

    def get_ext(self):
        filename = self.path_parse[-1]
        pos = filename.rfind(".")
        return filename[pos + 1:] if pos != -1 else ""

    def get_file_type(self):
        return self.file_type

    def get_sequence_delimiter(self):
        return self.read_delimiter


class PairFiles:
    def __init__(self):
        self.reset()

    def reset(self):
        self.pair_type = PairType.UnknownPair
        self.identify = ""
        self.open_correct = False
        self.first = SingleEndFile()
        self.second = SingleEndFile()

    def init(self, path_end1, path_end2):
        self.reset()

        init1 = self.first.init(path_end1)
        if init1:
            self.pair_type = PairType.SingleEnd
            init2 = self.second.init(path_end2)
            if init2:
                self.pair_type = PairType.PairedEnd
        else:
            init1 = self.first.init(path_end2)
            if init1:
                self.pair_type = PairType.SingleEnd
            else:
                self.pair_type = PairType.UnknownPair

        if init1:
            if self.pair_type == PairType.SingleEnd:
                self.identify = self.first.get_filename()
            elif self.pair_type == PairType.PairedEnd:
                self.identify = self.first.get_filename_without_ext()

        self.open_correct = init1
        if path_end1 and path_end2 and self.pair_type == PairType.PairedEnd:
            if self.first.get_file_type() != self.second.get_file_type():
                print("Different file type Fasta and Fastq")
                self.open_correct = False

        return self.open_correct

    def is_correct(self):
        return self.open_correct

    def get_pair_type(self):
        return self.pair_type

    def get_file_type(self):
        return self.first.get_file_type()

    def get_identify(self):
        return self.identify


# # Example usage:
# single_end_file = SingleEndFile()
# single_end_file.init("path/to/file.fasta")
# print(single_end_file.get_filename())

# pair_files = PairFiles()
# pair_files.init("path/to/file1.fasta", "path/to/file2.fasta")
# print(pair_files.get_identify())
