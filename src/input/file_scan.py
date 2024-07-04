class FileScan:
    def __init__(self):
        pass

    def scan(self, filename):
        lines = []
        try:
            with open(filename, 'r') as file:
                lines = file.readlines()
            return [line.strip() for line in lines]
        except IOError:
            return []
