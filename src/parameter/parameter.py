class Parameter:
    def __init__(self):
        self.spaced_qmers = []

    def init(self, param_file):
        try:
            with open(param_file, 'r') as file:
                # Additional initialization logic...
                pass
            return True
        except IOError:
            print(f"Unable to open parameter file: {param_file}")
            return False

    def add_spaced_qmer(self, qmer1, qmer2):
        self.spaced_qmers.append((qmer1, qmer2))
