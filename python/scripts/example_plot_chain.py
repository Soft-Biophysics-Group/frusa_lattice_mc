import config as cfg
import plotting.plot_chain as pc
import matplotlib.pyplot as plt
import numpy as np

_, orientations = np.loadtxt(cfg.data_path/'structures/final_structure.dat',
                             dtype=int)
fig, ax = plt.subplots()
pc.plot_chain(orientations, ax)
plt.show()
