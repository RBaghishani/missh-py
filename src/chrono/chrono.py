import time

class Chrono:
    def __init__(self):
        self.start_time = None
        self.end_time = None

    def start(self):
        self.start_time = time.time()

    def stop(self):
        self.end_time = time.time()

    def elapsed(self):
        if self.start_time is not None and self.end_time is not None:
            return self.end_time - self.start_time
        return 0.0

    def elapsed_str(self):
        elapsed_seconds = self.elapsed()
        hours = int(elapsed_seconds // 3600)
        elapsed_seconds -= hours * 3600
        minutes = int(elapsed_seconds // 60)
        elapsed_seconds -= minutes * 60
        seconds = int(elapsed_seconds)
        return f"{hours:02}:{minutes:02}:{seconds:02}"
