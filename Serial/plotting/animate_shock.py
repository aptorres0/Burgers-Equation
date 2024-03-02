#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation

# input data parameters
x_final = 2.0
#T_final = 0.30
T_final = 3.0

Z = np.loadtxt("/Users/ajpt/School/CompPhys-7411/Project-1/build/laxwendroff-c-output.txt", float)  # reads matrix u_ij
# print(Z.shape)
# print(Z)

fig, ax = plt.subplots(figsize=(5, 3))
ax.set(xlim=(0, x_final), ylim=(-5, 5))
ax.set_xlabel("x")
ax.set_ylabel("u(x,t)")

x = np.linspace(0, x_final, Z.shape[1])
t = np.linspace(0, T_final, Z.shape[0])
X2, T2 = np.meshgrid(x, t)

line = ax.plot(x, Z[0, :], color="k", lw=2)[0]


def animate(i):
    line.set_ydata(Z[i, :])


anim = FuncAnimation(fig, animate, interval=50, frames=len(t) - 1, repeat=False)

plt.draw()
plt.show()
