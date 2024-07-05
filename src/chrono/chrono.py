import time

class Chrono:
    def __init__(self):
        pass

    def exe(self, obj, function, param, ret_container):
        start_time = time.time()
        ret_container[0] = getattr(obj, function)(param)
        end_time = time.time()
        return (start_time, end_time)

    def exe_no_return(self, obj, function, param):
        start_time = time.time()
        getattr(obj, function)(param)
        end_time = time.time()
        return (start_time, end_time)

    def exe_two_params(self, obj, function, param1, param2):
        start_time = time.time()
        getattr(obj, function)(param1, param2)
        end_time = time.time()
        return (start_time, end_time)
