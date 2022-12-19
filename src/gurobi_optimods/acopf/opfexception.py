class OPFException(Exception):
    "Gurobi exception class"

    def _get_message(self):
        return self._message

    def _set_message(self, message):
        self._message = message

    message = property(_get_message, _set_message)

    def __init__(self, message):
        self.message = message

    def __str__(self):
        return self.message