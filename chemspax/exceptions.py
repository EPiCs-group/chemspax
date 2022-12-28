# -*- coding: utf-8 -*-
#                                                     #
#  __authors__ = Adarsh Kalikadien & Vivek Sinha      #
#  __institution__ = TU Delft                         #
#                                                     #


class RotationMatrixError(Exception):
    """Exception raised when c = -1.0
    Checks if angle with normal-vector and v_x != 0 or 180 else c = -1.0 which gives 1/0 in rotation matrix.
    """

    def __init__(self):
        super().__init__(
            "c = -1, which gives 1/0=undefined in rotation matrix. Exiting program"
        )


class InvalidRecursiveOrInitialArgumentError(Exception):
    def __init__(self):
        super().__init__(
            "recursive_or_initial argument is invalid. Exiting program"
        )
