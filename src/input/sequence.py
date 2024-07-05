class Sequence:
    def __init__(self):
        self.index_file = 0
        self.header = ""
        self.sequence = ""
        self.header_quality = ""
        self.quality = ""
        self.id = ""
        self.flag_end = ""

    def get_index_file(self):
        return self.index_file

    def set_index_file(self, index_file):
        self.index_file = index_file

    def get_header(self):
        return self.header

    def append_header(self, header, parser=None):
        self.header += header
        if parser:
            parser()

    def get_sequence(self):
        return self.sequence

    def append_sequence(self, sequence):
        self.sequence += sequence

    def get_header_quality(self):
        return self.header_quality

    def append_header_quality(self, header_quality):
        self.header_quality += header_quality

    def get_quality(self):
        return self.quality

    def append_quality(self, quality):
        self.quality += quality

    def get_id(self):
        return self.id

    def set_id(self, id):
        self.id = id

    def get_flag_end(self):
        return self.flag_end

    def set_flag_end(self, flag_end):
        self.flag_end = flag_end

    def is_sequence_all_n(self):
        return all(c == 'N' for c in self.sequence)

    def have_sequence_percent_n(self, perc):
        cont_n = self.sequence.count('N')
        perc_n = cont_n / len(self.sequence) if self.sequence else 0
        return perc_n >= perc

    def parser1(self):
        copy_header = self.header
        delimiter = ">r"
        pos = copy_header.find(delimiter)
        copy_header = copy_header[pos + len(delimiter):]
        delimiter = "."
        pos = copy_header.find(delimiter)
        
        self.set_id(copy_header[:pos])
        copy_header = copy_header[pos + len(delimiter):]
        self.set_flag_end(copy_header[0])

    def parser2(self):
        copy_header = self.header
        delimiter = "@"
        pos = copy_header.find(delimiter)
        copy_header = copy_header[pos + len(delimiter):]
        delimiter = "_"
        pos = copy_header.find(delimiter)
        
        self.set_id(copy_header[:pos])
        copy_header = copy_header[pos + len(delimiter):]
        self.set_flag_end(str(len(copy_header) - 1))