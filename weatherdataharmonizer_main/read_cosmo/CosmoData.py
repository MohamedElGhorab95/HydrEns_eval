#   @author: Michael Wagner
#   @organization: TU Dresden
#   @contact: michael.wagner@tu-dresden.de


class CosmoData:
    """
    CosmoData class contains all relevant data and metadata of Cosmo formatted data.
    """
    def __init__(self,
                 data=None,
                 metadata=None):
        """
        Initialize CosmoData class.

        :param data: 1D array of data
        :type data: numpy.ndarray
        :param metadata: describing metadata from icon data
        :type metadata: read_cosmo.Metadata.Metadata
        """
        self.data = data
        self.metadata = metadata
