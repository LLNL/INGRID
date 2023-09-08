
class RegionEntered(Exception):
    """
    Exception for identifying when we enter topological regions
    """
    def __init__(self, message, region):
        self.message = message
        self.region = region

class TkInitializationSuccess(Exception):
    """
    Exception for identifying we have reached the end of IngridGUI.__init__
    """
    pass