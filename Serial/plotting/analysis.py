#!/usr/bin/env python3

from math import ceil

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib.animation import FuncAnimation


def plot_check3():
    # Setup figure for plotting
    fig = plt.figure("breakdown_final_beta", constrained_layout=True)

    # Start out by plotting both methods using beta = 0.1 for tmax = 0.15
    Z_lf = np.loadtxt("../results/lf-b0.1-t0.3.txt", float)
    Z_lw = np.loadtxt("../results/lf-b1.2-t0.3.txt", float)
    num_points = Z_lf.shape[0]

    xmax = 2.0  # x [unitless]
    tmax = 0.3  # time [unitless]
    t = np.linspace(0, tmax, Z_lf.shape[0])
    x = np.linspace(0, xmax, Z_lf.shape[1])
    X, T = np.meshgrid(x, t)

    # create 2x1 subfigs
    subfigs = fig.subfigures(nrows=2, ncols=1)
    for row, subfig in enumerate(subfigs):
        if row == 0:
            subfig.suptitle(r"Leapfrog method for $\beta = 0.1$")
        if row == 1:
            subfig.suptitle(r"Leapfrog method for $\beta = 1.2$")

        # create 1x2 subplots per subfig
        axs = subfig.subplots(nrows=1, ncols=2)
        for col, ax in enumerate(axs):
            ax.grid()
            if row == 0 and col == 0:
                ax.set_xlabel("x")
                ax.set_ylabel(r"u(x,t)")
                ax.set_xlim([0, 2])

                ax.plot(
                    x,
                    Z_lf[0, :],
                    label=f"t = {0}",
                    antialiased=False,
                )

                point_n = 40
                t_point = (tmax / 150) * point_n
                ax.plot(
                    x,
                    Z_lf[point_n, :],
                    label=f"t = {t_point:0.3f}",
                    antialiased=False,
                )

                point_n = 57
                t_point = (tmax / 150) * point_n
                ax.plot(
                    x,
                    Z_lf[point_n, :],
                    label=f"t = {t_point:0.3f}",
                    antialiased=False,
                )

                point_n = 65
                t_point = (tmax / 150) * point_n
                ax.plot(
                    x,
                    Z_lf[point_n, :],
                    label=f"t = {t_point:0.3f}",
                    antialiased=False,
                )

                ax.legend(loc="upper right")
            if row == 1 and col == 0:
                ax.set_xlabel("x")
                ax.set_ylabel(r"u(x,t)")
                ax.set_xlim([0, 2])

                ax.plot(
                    x,
                    Z_lw[0, :],
                    label=f"t = {0}",
                    antialiased=False,
                )

                point_n = 2
                t_point = (tmax / 12) * point_n
                ax.plot(
                    x,
                    Z_lw[point_n, :],
                    label=f"t = {t_point:0.3f}",
                    antialiased=False,
                )

                point_n = 4
                t_point = (tmax / 12) * point_n
                ax.plot(
                    x,
                    Z_lw[point_n, :],
                    label=f"t = {t_point:0.3f}",
                    antialiased=False,
                )

                point_n = 5
                t_point = (tmax / 12) * point_n
                ax.plot(
                    x,
                    Z_lw[point_n, :],
                    label=f"t = {t_point:0.3f}",
                    antialiased=False,
                )

                point_n = 6
                t_point = (tmax / 12) * point_n
                ax.plot(
                    x,
                    Z_lw[point_n, :],
                    label=f"t = {t_point:0.3f}",
                    antialiased=False,
                )

                point_n = 7
                t_point = (tmax / 12) * point_n
                ax.plot(
                    x,
                    Z_lw[point_n, :],
                    label=f"t = {t_point:0.3f}",
                    antialiased=False,
                )

                ax.legend(loc="upper right")
            if row == 0 and col == 1:
                ax.set_xlabel("x")
                ax.set_xlim([0, 2])

                point_n = 60
                t_point = (tmax / 150) * point_n
                ax.plot(
                    x,
                    Z_lf[point_n, :],
                    label=f"t = {t_point:0.3f}",
                    antialiased=False,
                )

                point_n = 70
                t_point = (tmax / 150) * point_n
                ax.plot(
                    x,
                    Z_lf[point_n, :],
                    label=f"t = {t_point:0.3f}",
                    antialiased=False,
                )

                point_n = 80
                t_point = (tmax / 150) * point_n
                ax.plot(
                    x,
                    Z_lf[point_n, :],
                    label=f"t = {t_point:0.3f}",
                    antialiased=False,
                    alpha=0.5,
                )

                ax.legend(loc="upper right")

            if row == 1 and col == 1:
                ax.set_xlabel("x")
                ax.set_xlim([0, 2])

                point_n = 8
                t_point = (tmax / 12) * point_n
                ax.plot(
                    x,
                    Z_lw[point_n, :],
                    label=f"t = {t_point:0.3f}",
                    antialiased=False,
                )

                point_n = 9
                t_point = (tmax / 12) * point_n
                ax.plot(
                    x,
                    Z_lw[point_n, :],
                    label=f"t = {t_point:0.3f}",
                    antialiased=False,
                )

                point_n = 10
                t_point = (tmax / 12) * point_n
                ax.plot(
                    x,
                    Z_lw[point_n, :],
                    label=f"t = {t_point:0.3f}",
                    antialiased=False,
                    alpha=0.5,
                )

                ax.legend(loc="upper right")


plot_check3()

plt.show()
