#   @author: Michael Wagner
#   @organization: TU Dresden
#   @contact: michael.wagner@tu-dresden.de
from met_entities.VariableDescription import DataDescription


class LonLatTime:
    """
    LonLatTime class is used as a container for coordinate- or time-related data to facilitate the export to netCDF.
    """
    def __init__(self,
                 data=None,
                 data_description=None):
        """
        Initialize LonLatTime object.

        :param data: values as a n-dimensional matrix
        :type data: numpy.ndarray, optional
        :param data_description: metadata for the data
        :type data_description: VariableDescription.DataDescription or
            VariableDescription.TimeDescription, optional
        """
        self.data = data
        self.data_description = data_description

    def __copy__(self):
        """
        Provide deep copy of LonLatTime instance.

        :return: copied instance of this LonLatTime object
        :rtype: LonLatTime
        """
        copy_instance = LonLatTime(data=self.data)
        copy_instance.data_description = DataDescription(
            fill_value=self.data_description.fill_value,
            units=self.data_description.units,
            scale_factor=self.data_description.scale_factor,
            coordinate_system=self.data_description.coordinate_system,
            long_name=self.data_description.long_name,
            standard_name=self.data_description.standard_name,
            time_note=self.data_description.time_note,
            eps_note=self.data_description.eps_note)

        return copy_instance
