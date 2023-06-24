#   @author: Michael Wagner
#   @organization: TU Dresden
#   @contact: michael.wagner@tu-dresden.de


class IconData:
    """
    IconData class contains all relevant data and metadata of Icon formatted data.
    """
    def __init__(self,
                 data=None,
                 metadata=None):
        """
        Initialize IconData class.

        :param data: 1D array of data
        :type data: numpy.ndarray
        :param metadata: describing metadata from icon data
        :type metadata: read_icon.Metadata.Metadata
        """
        self.data = data
        self.metadata = metadata
