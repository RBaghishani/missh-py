import os

class PositionRead:
    def __init__(self, start=None, end=None):
        self.start = start
        self.end = end

class FileScan:
    def __init__(self):
        self.file = SingleEndFile()
        self.stream = None
        self.pos_read = []

    def __del__(self):
        if self.stream:
            self.stream.close()

    def __copy__(self):
        new_copy = FileScan()
        new_copy.init(self.file.get_path())
        return new_copy

    def init(self, path):
        ret = self.file.init(path)
        if ret:
            ret &= self.get_stream_pos_read(self.file.get_path())
        return ret

    def get_sequence_number(self):
        return len(self.pos_read)

    def get_sequence_with_index(self, index, read, parser_header=None):
        if index >= len(self.pos_read):
            return False
        if not self.stream:
            self.stream = open(self.file.get_path(), 'r')

        actual = self.stream.tell()
        if actual != self.pos_read[index].start:
            self.stream.seek(self.pos_read[index].start)

        parse_header = parser_header is not None
        is_line_quality = False

        while actual != self.pos_read[index].end:
            line = self.stream.readline().strip()
            if not line:
                break

            if line[0] == self.file.get_sequence_delimiter():
                if not is_line_quality:
                    read.set_index_file(index)
                    if parse_header:
                        read.append_header(line, parser_header)
                    else:
                        read.append_header(line)
                else:
                    read.append_quality(line)
                    is_line_quality = False
            else:
                if self.file.get_file_type() == FileType.Fastq and line[0] == '+':
                    read.append_header_quality(line)
                    is_line_quality = True
                else:
                    if not is_line_quality:
                        read.append_sequence(line)
                    else:
                        read.append_quality(line)
                        is_line_quality = False
            actual = self.stream.tell()
        return True

    def get_stream_pos_read(self, path):
        with open(path, 'r') as file:
            is_line_quality = False
            before_call_get_line = file.tell()

            while True:
                line = file.readline()
                if not line:
                    break

                if line[0] == self.file.get_sequence_delimiter():
                    if not is_line_quality:
                        if self.pos_read:
                            self.pos_read[-1].end = before_call_get_line
                        pos = PositionRead(start=before_call_get_line)
                        self.pos_read.append(pos)
                    else:
                        is_line_quality = False
                else:
                    if self.file.get_file_type() == FileType.Fastq and line[0] == '+':
                        is_line_quality = True
                    else:
                        is_line_quality = False

                before_call_get_line = file.tell()

            if self.pos_read:
                self.pos_read[-1].end = before_call_get_line

        self.stream = open(path, 'r')
        return True

    def get_file(self):
        return self.file

    def is_correct(self):
        return self.file.is_correct()

    def reset(self):
        self.file.reset()
        if self.stream:
            self.stream.close()
        self.stream = None
        self.pos_read = []


class FilesScan:
    def __init__(self):
        self.reset()

    def __del__(self):
        pass

    def __copy__(self):
        new_copy = FilesScan()
        new_copy.init(self.first.get_file().get_path(), self.second.get_file().get_path())
        return new_copy

    def reset(self):
        self.pair_type = PairType.UnknownPair
        self.identify = ""
        self.open_correct = False
        self.first = FileScan()
        self.second = FileScan()

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
                self.identify = self.first.get_file().get_filename()
            elif self.pair_type == PairType.PairedEnd:
                self.identify = self.first.get_file().get_filename_without_ext()

        self.open_correct = init1
        if path_end1 and path_end2 and self.pair_type == PairType.PairedEnd:
            if self.first.get_sequence_number() != self.second.get_sequence_number():
                self.open_correct = False
                print(f"Different number of read in file {self.first.get_file().get_path()} and {self.second.get_file().get_path()}")
            if self.first.get_file().get_file_type() != self.second.get_file().get_file_type():
                self.open_correct = False
                print("Different file type Fasta and Fastq")
        return self.open_correct

    def get_paired_reads_number(self):
        return self.first.get_sequence_number()

    def get_sequences_number(self):
        return self.first.get_sequence_number() + self.second.get_sequence_number()

    def get_identify(self):
        return self.identify

    def get_pair_type(self):
        return self.pair_type

    def get_file_type(self):
        return self.first.get_file().get_file_type()

    def is_correct(self):
        return self.open_correct
