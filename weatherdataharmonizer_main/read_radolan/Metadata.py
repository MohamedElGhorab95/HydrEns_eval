#   @author: Michael Wagner
#   @organization: TU Dresden
#   @contact: michael.wagner@tu-dresden.de


class Metadata:
    """
    Metadata class contains relevant metadata of radolan formatted data.
    """
    def __init__(self,
                 product=None,
                 datum=None,
                 datum_iso=None,
                 product_length=None,
                 format_version=None,
                 version=None,
                 interval_length=None,
                 ncols=None,
                 nrows=None,
                 coord_ll=None,
                 prediction_time=None,
                 module_flags=None,
                 quantification_kind=None,
                 sites=None):
        """
        Initialize Metadata class.

        :param product: abbreviation of the product
        :type product: str
        :param datum: UTC datum/time of observation (yy mm dd HH MM)
        :type datum: list
        :param datum_iso: ISO 8601 compliant datetime of observation
        :type datum_iso: str
        :param product_length: product length in byte
        :type product_length: str
        :param format_version: format version; 0 for mixed version (100 and 128 km), 1 for 100 km, 2 for 128 km,
            3 for 150 km
        :type format_version: int
        :param version: software version
        :type version: str
        :param interval_length: interval length in minutes
        :type interval_length: int
        :param ncols: number of columns
        :type ncols: int
        :param nrows: number of rows
        :type nrows: int
        :param coord_ll: coordinates of lower left corner
        :type coord_ll: list
        :param prediction_time: prediction time in minutes after observation
        :type prediction_time: int
        :param module_flags: module flags
        :type module_flags: int
        :param quantification_kind: defined quantification
        :type quantification_kind: int
        :param sites: abbreviated radar sites used for the merged product
        :type sites: str
        """
        if datum is None:
            datum = []
        self.product = product
        self.datum = datum
        self.datum_iso = datum_iso
        self.product_length = product_length
        self.format_version = format_version
        self.version = version
        self.interval_length = interval_length
        self.ncols = ncols
        self.nrows = nrows
        self.coord_ll = coord_ll
        self.prediction_time = prediction_time
        self.module_flags = module_flags
        self.quantification_kind = quantification_kind
        self.sites = sites
        