"""
A wrapper for the lisem agnostic implementation of the dynamic Manning's value fpr partially submerged vegetation

Usage:
calcManning.functiontype = 1  # 1: Oberle function
NN = calcManning(NN=0.035, Whr=0.1, PH=0.2)

"""

import ctypes
from pathlib import Path
import os

class JLUManning:
    """
    Wrapper class
    """
    def compile(self):
        """
        Compiles the code if the .so library does not exist
        """
        libname = self.cdll_file
        srcname = self.home / 'jlu_manning.cpp'
        os.system(f'g++ -shared -o {libname} {srcname}')

    def __init__(self, functiontype: int):
        self.functiontype = functiontype
        self.home = Path(__file__).parent.absolute()
        self.cdll_file = self.home / 'libjlu_manning.so'
        if not self.cdll_file.exists():
            self.compile()
        self.cdll = ctypes.CDLL(self.cdll_file)
        self.cdll.calcManningDispatch.restype = ctypes.c_double

    def __call__(self, NN: float, WHr: float=0, PH:float=0, coverc: float=0):
        """
        Calls the wrpped C++ function
        :param NN: Manning's N : scope of it depends on the function
        :param WHr: Water Height (in m)
        :param PH: Plant height (in m)
        :param coverc: Plant coverage
        :return:
        """
        return self.cdll.calcManningDispatch(
            ctypes.c_int(self.functiontype),
            ctypes.c_double(NN),
            ctypes.c_double(WHr),
            ctypes.c_double(PH),
            ctypes.c_double(coverc)
        )


calcManning = JLUManning(1)
