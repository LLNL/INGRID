from ingrid.core.geometry.boundaries import ClosedBoundary

class Patch:
    boundary: ClosedBoundary
    def __init__(self, boundary: ClosedBoundary) -> None:
        self.boundary = boundary

    