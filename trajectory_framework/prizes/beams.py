import numpy as np
import mpl_toolkits.mplot3d
import matplotlib.pyplot as pp
import math
import time

from utils.options import get_options


class get_beams():

    @staticmethod
    def full_path(opts):
        """create a set of preliminary beams based on the original 
        set of nodes
        """
        num_pts = opts.num_beams * 2
        indices = np.arange(0, num_pts, dtype=float) + 0.5

        phi = np.arccos(1 - 2*indices/num_pts)
        theta = np.pi * (1 + 5**0.5) * indices

        # update phi and theta
        points = list(zip(phi, theta))
        output = [i for i in points if i[0] < np.pi/2]
        phi, theta = list(zip(*output))

        print(max(phi))

        return phi, theta

    @staticmethod
    def fixed_path(opts):
        """The fixed path is extracted the full body path that is
        currently being used clinically at CHUM, which is a fixed 
        set of nodes and direction points towards the body
        """
        pass

    @staticmethod
    def plot_beams(phi, theta):
        x, y, z = np.cos(theta) * np.sin(phi), np.sin(theta) * \
            np.sin(phi), np.cos(phi)
        pp.figure().add_subplot(111, projection='3d').scatter(x, y, z)
        pp.show()


if __name__ == "__main__":
    run(get_options())
