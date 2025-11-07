#-----------------------------------------------------------------------
class Errors:
    """
    Class for defining error measures for optimiser. 

    Author: Z. Nikolaou
    """

    @staticmethod
    def msqe(a, b): return np.square(np.subtract(a, b)).mean()

    @staticmethod
    def mase(a, b): return np.abs(np.subtract(a, b)).mean()
#-----------------------------------------------------------------------

