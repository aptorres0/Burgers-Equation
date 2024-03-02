#!/usr/bin/env python3

from math import ceil

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib.animation import FuncAnimation

# Setup plotting
# plt.style.use("seaborn-whitegrid")
# Set Matplotlib defaults
# plt.rc("figure", autolayout=True)
# plt.rc(
#    "axes",
#    labelweight="bold",
#    labelsize="large",
#    titleweight="bold",
#    titlesize=18,
#    titlepad=10,
# )


def plot_b():
    # Setup figure for plotting
    fig = plt.figure("plot_b", constrained_layout=True)

    # Start out by plotting both methods using beta = 0.1 for tmax = 0.15
    data_lf = np.loadtxt("../results/lf-b0.1-t0.15.txt", float)
    data_lw = np.loadtxt("../results/lw-b0.1-t0.15.txt", float)

    xmax = 2.0  # x [unitless]
    tmax = 0.15  # time [unitless]
    t = np.linspace(0, tmax, data_lf.shape[0])
    x = np.linspace(0, xmax, data_lf.shape[1])
    X, T = np.meshgrid(x, t)

    fig2 = plt.figure("wire", constrained_layout=True)
    ax2 = fig2.add_subplot(111, projection="3d")
    ax2.plot_wireframe(X, T, data_lw)

    ax = fig.add_subplot(231, projection="3d")
    ax.plot_surface(X, T, data_lf, cmap=cm.Spectral, antialiased=False)
    zmin, _ = ax.get_zlim()
    # ax.contour(X, T, data_lf, zdir="z", offset=zmin)
    ax.set(
        title=r"Leapfrog Method for $\beta = 0.1$",
        xlim=[0, xmax],
        ylim=[0, tmax],
        xlabel="space",
        ylabel="time",
        zlabel=r"u(x,t)",
    )

    ax = fig.add_subplot(232, projection="3d")
    ax.plot_surface(X, T, data_lw, cmap=cm.Spectral, lw=0, antialiased=False)
    ax.set(
        title=r"Lax-Wendroff Method for $\beta = 0.1$",
        xlim=[0, xmax],
        ylim=[0, tmax],
        xlabel="space",
        ylabel="time",
    )

    # Now plot both methods for beta = 0.1 and tmax = 0.3
    data_lf = np.loadtxt("../results/lf-b0.1-t0.3.txt", float)
    data_lw = np.loadtxt("../results/lw-b0.1-t0.3.txt", float)

    xmax = 2.0  # x [unitless]
    tmax = 0.3  # time [unitless]
    t = np.linspace(0, tmax, data_lf.shape[0])
    x = np.linspace(0, xmax, data_lf.shape[1])
    X, T = np.meshgrid(x, t)

    ax = fig.add_subplot(233, projection="3d")
    ax.plot_surface(X, T, data_lf, cmap=cm.Spectral, lw=0, antialiased=False)
    ax.set(
        title=r"Leapfrog Method for $\beta = 0.1$",
        xlim=[0, xmax],
        ylim=[0, tmax],
        xlabel="space",
        ylabel="time",
    )

    ax = fig.add_subplot(234, projection="3d")
    ax.plot_surface(X, T, data_lw, cmap=cm.Spectral, lw=0, antialiased=False)
    ax.set(
        title=r"Lax-Wendroff for $\beta = 0.1$",
        xlim=[0, xmax],
        ylim=[0, tmax],
        xlabel="space",
        ylabel="time",
    )

    # Now plot both methods for beta = 0.1 and tmax = 0.2
    data_lf = np.loadtxt("../results/lf-b0.1-t0.2.txt", float)
    data_lw = np.loadtxt("../results/lw-b0.1-t0.2.txt", float)

    xmax = 2.0  # x [unitless]
    tmax = 0.2  # time [unitless]
    t = np.linspace(0, tmax, data_lf.shape[0])
    x = np.linspace(0, xmax, data_lf.shape[1])
    X, T = np.meshgrid(x, t)

    ax = fig.add_subplot(235, projection="3d")
    ax.plot_surface(X, T, data_lf, cmap=cm.Spectral, lw=0, antialiased=False)
    ax.set(
        title=r"Leapfrog Method for $\beta = 0.1$",
        xlim=[0, xmax],
        ylim=[0, tmax],
        xlabel="space",
        ylabel="time",
    )

    ax = fig.add_subplot(236, projection="3d")
    ax.plot_surface(X, T, data_lw, cmap=cm.Spectral, lw=0, antialiased=False)
    ax.set(
        title=r"Lax-Wendroff for $\beta = 0.1$",
        xlim=[0, xmax],
        ylim=[0, tmax],
        xlabel="space",
        ylabel="time",
    )


def plot_b2():

    xmax = 2.0  # x [unitless]

    fig = plt.figure("lf_b01_t25_wire", constrained_layout=True)
    data_lf = np.loadtxt("../results/lf-b0.1-t0.3.txt", float)
    data_lw = np.loadtxt("../results/lw-b0.1-t0.3.txt", float)
    tmax = 0.3  # time [unitless]

    n_point = int((0.25 / tmax) * len(data_lf))

    t = np.linspace(0, 0.25, n_point)
    # t = np.linspace(0, tmax, data_lw.shape[0])
    x = np.linspace(0, xmax, data_lw.shape[1])

    X, T = np.meshgrid(x, t)

    ax = fig.add_subplot(111, projection="3d")
    surf = ax.plot_wireframe(X, T, data_lf[:n_point, :], antialiased=False)
    ax.set(
        title=r"Leapfrog Method for $\beta = 0.1$",
        xlim=[0, xmax],
        ylim=[0, 0.25],
        xlabel="space",
        ylabel="time",
        zlabel=r"u(x,t)",
    )

    # zmin, _ = ax.get_zlim()
    # ax.contour(X, T, data_lf, zdir="z", offset=zmin)
    ax.set_ylabel("time", labelpad=10)
    ax.set_xlabel("space", labelpad=10)
    ax.set_zlabel(r"u(x,t)", labelpad=10)

    print(len(data_lf))


def plot_check():
    # Setup figure for plotting
    fig = plt.figure(figsize=(8, 8))

    # Start out by plotting both methods using beta = 0.1 for tmax = 0.15
    Z_lf = np.loadtxt("../results/lf-b0.1-t0.3.txt", float)
    Z_lw = np.loadtxt("../results/lw-b0.1-t0.3.txt", float)
    num_points = Z_lf.shape[0]

    xmax = 2.0  # x [unitless]
    tmax = 0.3  # time [unitless]
    t = np.linspace(0, tmax, Z_lf.shape[0])
    x = np.linspace(0, xmax, Z_lf.shape[1])
    X, T = np.meshgrid(x, t)

    ax_lf = fig.add_subplot(211)
    ax_lf.set_title(r"Leapfrog method for $\beta = 0.1$")
    ax_lf.set_xlabel("x")
    ax_lf.set_ylabel("u(x,t)")

    ax_lw = fig.add_subplot(212)
    ax_lw.set_title(r"Lax-Wendroff method for $\beta = 0.1$")
    ax_lw.set_xlabel("x")
    ax_lw.set_ylabel("u(x,t)")

    # Plot at t = 0
    ax_lf.plot(x, Z_lf[0, :], label=f"t = {0}")
    ax_lw.plot(x, Z_lw[0, :], label=f"t = {0}")

    # Plot at quarterway point
    point_n = num_points // 4
    t_point = tmax / 4
    ax_lf.plot(x, Z_lf[point_n, :], label=f"t = {t_point:0.3f}")
    ax_lw.plot(x, Z_lw[point_n, :], label=f"t = {t_point:0.3f}")

    # Plot at halfway point
    point_n = num_points // 2
    t_point = tmax / 2
    ax_lf.plot(x, Z_lf[point_n, :], label=f"t = {t_point:0.3f}")
    ax_lw.plot(x, Z_lw[point_n, :], label=f"t = {t_point:0.3f}")

    # Plot at 3/4 tmax
    point_n = 3 * ceil(num_points / 4)
    t_point = tmax * 3 / 4
    # ax_lf.plot(x, Z_lf[point_n, :], label=f"t = {t_point:0.3f}")
    ax_lw.plot(x, Z_lw[point_n, :], label=f"t = {t_point:0.3f}")

    # Plot at t = tmax
    # ax_lf.plot(x, Z_lf[num_points - 1, :], label=f"t = {tmax:0.3f}")
    ax_lw.plot(x, Z_lw[num_points - 1, :], label=f"t = {tmax:0.3f}")

    ax_lf.legend(loc="upper right")
    ax_lw.legend(loc="upper right")


def plot_check2():
    # Setup figure for plotting
    fig = plt.figure("breakdown", constrained_layout=True)

    # Start out by plotting both methods using beta = 0.1 for tmax = 0.15
    Z_lf = np.loadtxt("../results/lf-b0.1-t0.3.txt", float)
    Z_lw = np.loadtxt("../results/lw-b0.1-t0.3.txt", float)
    num_points = Z_lf.shape[0]

    xmax = 2.0  # x [unitless]
    tmax = 0.3  # time [unitless]
    t = np.linspace(0, tmax, Z_lf.shape[0])
    x = np.linspace(0, xmax, Z_lf.shape[1])
    X, T = np.meshgrid(x, t)

    ax_lf = fig.add_subplot(221)
    ax_lf.set_title(r"Leapfrog method for $\beta = 0.1$")
    ax_lf.set_xlabel("x")
    ax_lf.set_ylabel("u(x,t)")

    ax_lw = fig.add_subplot(223)
    ax_lw.set_title(r"Lax-Wendroff method for $\beta = 0.1$")
    ax_lw.set_xlabel("x")
    ax_lw.set_ylabel("u(x,t)")

    ax_lf2 = fig.add_subplot(222)
    ax_lf2.set_title(r"Leapfrog method for $\beta = 0.1$")
    ax_lf2.set_xlabel("x")
    ax_lf2.set_ylabel("u(x,t)")

    ax_lw2 = fig.add_subplot(224)
    ax_lw2.set_title(r"Lax-Wendroff method for $\beta = 0.1$")
    ax_lw2.set_xlabel("x")
    ax_lw2.set_ylabel("u(x,t)")

    # Plot at t = 0
    ax_lf.plot(x, Z_lf[0, :], label=f"t = {0}")
    ax_lw.plot(x, Z_lw[0, :], label=f"t = {0}")

    point_n = 40
    t_point = (tmax / 150) * point_n
    ax_lf.plot(x, Z_lf[point_n, :], label=f"t = {t_point:0.3f}")
    ax_lw.plot(x, Z_lw[point_n, :], label=f"t = {t_point:0.3f}")

    point_n = 59
    t_point = (tmax / 150) * point_n
    ax_lw.plot(x, Z_lw[point_n, :], label=f"t = {t_point:0.3f}")

    point_n = 60
    t_point = (tmax / 150) * point_n
    ax_lf2.plot(x, Z_lf[point_n, :], label=f"t = {t_point:0.3f}")
    ax_lw2.plot(x, Z_lw[point_n, :], label=f"t = {t_point:0.3f}")

    point_n = 70
    t_point = (tmax / 150) * point_n
    ax_lf2.plot(x, Z_lf[point_n, :], label=f"t = {t_point:0.3f}")
    ax_lw2.plot(x, Z_lw[point_n, :], label=f"t = {t_point:0.3f}")

    point_n = 80
    t_point = (tmax / 150) * point_n
    ax_lf2.plot(x, Z_lf[point_n, :], label=f"t = {t_point:0.3f}")
    ax_lw2.plot(x, Z_lw[point_n, :], label=f"t = {t_point:0.3f}")

    ax_lf.legend(loc="upper right")
    ax_lw.legend(loc="upper right")
    ax_lf2.legend(loc="upper right")
    ax_lw2.legend(loc="upper right")


def test():
    # Setup figure for plotting
    fig = plt.figure(figsize=((5, 3)))

    # Start out by plotting both methods using beta = 0.1 for tmax = 0.15
    Z_lf = np.loadtxt("../results/lf-b0.1-t0.3.txt", float)
    Z_lw = np.loadtxt("../results/lw-b0.1-t0.3.txt", float)

    xmax = 2.0  # x [unitless]
    tmax = 0.3  # time [unitless]
    t = np.linspace(0, tmax, Z_lf.shape[0])
    x = np.linspace(0, xmax, Z_lf.shape[1])
    X, T = np.meshgrid(x, t)

    ax = fig.add_subplot(111)
    ax.set_title("test")
    ax.set_xlabel("x")
    ax.set_ylabel("u(x,t)")

    # Plot at t = 0
    ax.plot(x, Z_lf[0, :], label=f"t = {0}")

    # Plot at halfway point
    num_points = Z_lf.shape[0]
    point_n = num_points // 4
    t_point = tmax / 4
    ax.plot(x, Z_lf[point_n, :], label=f"t = {t_point:0.3f}")

    # Plot at t = 0.15
    point_n = num_points // 2
    t_point = tmax / 2
    ax.plot(x, Z_lf[point_n, :], label=f"t = {t_point:0.3f}")

    ax.legend()


def plot_c():
    fig = plt.figure("frames", constrained_layout=True)
    # fig.suptitle('Figure title')

    xmax = 2.0  # x [unitless]
    tmax = 0.3  # time [unitless]

    # create 3x1 subfigs
    subfigs = fig.subfigures(nrows=3, ncols=1)
    for row, subfig in enumerate(subfigs):
        if row == 0:
            subfig.suptitle(r"Leapfrog method for $\beta = 0.1$")
            Z = np.loadtxt("../results/lf-b0.1-t0.3.txt", float)

        if row == 1:
            subfig.suptitle(r"Lax-Wendroff method for $\beta = 0.1$")
            Z = np.loadtxt("../results/lw-b0.1-t0.3.txt", float)

        if row == 2:
            subfig.suptitle(r"Lax-Wendroff method for $\beta = 1.2$")
            Z = np.loadtxt("../results/lw-b1.2-t0.3.txt", float)

        num_points = Z.shape[0]
        t = np.linspace(0, tmax, Z.shape[0])
        x = np.linspace(0, xmax, Z.shape[1])
        X, T = np.meshgrid(x, t)

        # create 1x3 subplots per subfig
        axs = subfig.subplots(nrows=1, ncols=3)
        for col, ax in enumerate(axs):
            if row == 0:
                if col == 0:
                    point_n = num_points // 3
                    t_point = tmax / 3
                    ax.plot(x, Z[point_n, :], label=f"t = {t_point:0.3f}")
                    ax.set_title(f"t = {t_point:0.3f}")
                    ax.set_xlabel(r"$x$")
                    ax.grid()
                if col == 1:
                    # Plot at quarterway point
                    t_point = 0.125
                    ax.plot(x, Z[63, :], label=f"t = {t_point:0.3f}")
                    ax.set_title(f"t = {t_point:0.3f}")
                    ax.set_xlabel(r"$x$")
                    ax.grid()
                if col == 2:
                    t_point = 0.12
                    ax.plot(x, Z[55, :], label=f"t = {t_point:0.3f}")
                    ax.set_title(f"t = {t_point:0.3f}")
                    ax.set_xlabel(r"$x$")
                    ax.grid()
                    # point_n = num_points // 2
                    # t_point = tmax / 2
                    # ax.plot(x, Z[point_n, :], label=f"t = {t_point:0.3f}")
                    # ax.set_title(f"t = {t_point:0.3f}")
                    # ax.set_xlabel(r"$x$")
                    # ax.grid()

            if row == 1:
                if col == 0:
                    point_n = num_points // 3
                    t_point = tmax / 3
                    ax.plot(x, Z[point_n, :], label=f"t = {t_point:0.3f}")
                    ax.set_title(f"t = {t_point:0.3f}")
                    ax.set_xlabel(r"$x$")
                    ax.grid()
                if col == 1:
                    # Plot at 0.125
                    t_point = 0.11111111
                    ax.plot(x, Z[58, :], label=f"t = {t_point:0.3f}")
                    ax.set_title(f"t = {t_point:0.3f}")
                    ax.set_xlabel(r"$x$")
                    ax.grid()
                if col == 2:
                    # Plot at halfway point
                    point_n = num_points // 2
                    t_point = tmax / 2
                    ax.plot(x, Z[point_n, :], label=f"t = {t_point:0.3f}")
                    ax.set_title(f"t = {t_point:0.3f}")
                    ax.set_xlabel(r"$x$")
                    ax.grid()

            if row == 2:
                if col == 0:
                    ax.plot(x, Z[0, :], label=f"t = {0}")
                    ax.plot()
                    ax.set_title(r"$t = 0$")
                    ax.set_xlabel(r"$x$")
                    ax.set_ylabel(r"$u(x,t)$")
                    ax.grid()
                if col == 1:
                    ax.plot(x, Z[0, :], label=f"t = {0}")
                    ax.plot()
                    ax.set_title(r"$t = 0$")
                    ax.set_xlabel(r"$x$")
                    ax.set_ylabel(r"$u(x,t)$")
                    ax.grid()
                    # Plot at quarterway point
                if col == 2:
                    point_n = num_points // 2
                    t_point = tmax / 2
                    ax.plot(x, Z[point_n, :], label=f"t = {t_point:0.3f}")
                    ax.set_title(f"t = {t_point:0.3f}")
                    ax.set_xlabel(r"$x$")
                    ax.grid()


def plot_c_ani_lf():
    # Setup figure for plotting
    # fig = plt.figure()
    # fig.figsize((5, 3))
    fig, ax = plt.subplots(figsize=(5, 3))

    # Start out by plotting both methods using beta = 0.1 for tmax = 0.15
    Z_lf = np.loadtxt("../results/lf-b0.1-t0.15.txt", float)
    Z_lw = np.loadtxt("../results/lw-b0.1-t0.15.txt", float)

    xmax = 2.0  # x [unitless]
    tmax = 0.15  # time [unitless]
    t = np.linspace(0, tmax, Z_lf.shape[0])
    x = np.linspace(0, xmax, Z_lf.shape[1])
    X, T = np.meshgrid(x, t)

    line_lf = ax.plot(x, Z_lf[0, :], color="k", lw=2)[0]

    def animate(i):
        line_lf.set_ydata(Z_lf[i, :])

    anim = FuncAnimation(fig, animate, interval=50, frames=len(t) - 1, repeat=False)

    return anim


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


plot_b2()
# plot_check3()

plt.show()
