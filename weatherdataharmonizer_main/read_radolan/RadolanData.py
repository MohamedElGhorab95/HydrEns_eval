#   @author: Michael Wagner
#   @organization: TU Dresden
#   @contact: michael.wagner@tu-dresden.de


class RadolanData:
    """
    RadolanData class contains all relevant data of radolan formatted data.
    """
    def __init__(self,
                 data=None,
                 metadata=None,
                 idx_clutter=None,
                 idx_nan=None):
        """
        Initialize RadolanData class.

        :param data: 2D matrix with radolan data
        :type data: numpy.ndarray
        :param metadata: describing metadata from radolan data
        :type metadata: read_radolan.Metadata.Metadata
        :param idx_clutter: 2D matrix with True/1 for clutter in radolan data und False/0 otherwise
        :type idx_clutter: numpy.ndarray
        :param idx_nan: 2D matrix with True/1 for nan in radolan data and False/0 otherwise
        :type idx_nan: numpy.ndarray
        """
        self.data = data
        self.metadata = metadata
        self.idx_clutter = idx_clutter
        self.idx_nan = idx_nan
