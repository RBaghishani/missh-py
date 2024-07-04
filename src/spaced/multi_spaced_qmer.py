class MultiSpacedQmer:
    def __init__(self):
        self.spaced_qmers = []

    def add_spaced_qmer(self, qmer):
        self.spaced_qmers.append(qmer)

    def get_spaced_qmers(self):
        return self.spaced_qmers
