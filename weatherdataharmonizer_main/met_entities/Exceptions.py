#   @author: Michael Wagner
#   @organization: TU Dresden
#   @contact: michael.wagner@tu-dresden.de


class RadarFileNotAvailable(Exception):
    """
    Raised if a radar observation file is not found.
    """
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


class ForecastFileNotAvailable(Exception):
    """
    Raised if already the first file of forecast dataset is not found.
    """
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)

class ForecastFileFormatError(Exception):
    """
    Raised if the format of a forecast file does not follow the rules.
    """
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)
