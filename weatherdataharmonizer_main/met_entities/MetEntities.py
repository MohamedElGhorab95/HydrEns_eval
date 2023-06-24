#   @author: Michael Wagner
#   @organization: TU Dresden
#   @contact: michael.wagner@tu-dresden.de

from abc import ABC, abstractmethod

from met_entities.LonLatTime import LonLatTime
from met_entities.GeoReferencedData import GeoReferencedData


class MetEntities(ABC):
    """
    MetEntities is an abstract class that defines methods for handling meteorological data.
    """

    def __init__(self,
                 time_value: LonLatTime = None,
                 forecast_value: LonLatTime = None,
                 gr_data=None,
                 eps_member: list = None,
                 short: str = False):
        self.time_value = time_value
        self.forecast_value = forecast_value
        self.gr_data = gr_data
        self.eps_member = eps_member
        self.short = short

    @abstractmethod
    def read_file(self, **kwargs):
        pass

    @abstractmethod
    def regrid(self, **kwargs):
        pass

    @abstractmethod
    def crop(self, **kwargs):
        pass

    @abstractmethod
    def export_netcdf(self, filename, **kwargs):
        pass

    @abstractmethod
    def export_netcdf_append(self, filename):
        pass
