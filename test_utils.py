import time

class Timer:
    def __init__(self, time_limit: float = 60):
        self._start_time = time.time()
        self._lap_time = time.time()
        self._time_limit = time_limit

    def time(self) -> float:
        return time.time() - self._start_time

    def lap_time(self) -> float:
        return time.time() - self._lap_time

    def start_lap(self):
        self._lap_time = time.time()

    def timed_out(self):
        return (time.time() - self._start_time) > self._time_limit