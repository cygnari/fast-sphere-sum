import numpy as np
import scipy as sp
import math
import matplotlib.pyplot as plt

Lx = 10
Ly = 10
Lz = 10

dx = 1
dy = 1
dz = 1

Nx = int(Lx / dx) + 1
Ny = int(Ly / dy) + 1
Nz = int(Lz / dz) + 1

K = 75
E = 100
lambd = 3 * K * (3 * K - E) / (9 * K - E)
mu = 3 * K * E / (9 * K - E)
sigma = lambd + 2 * mu
eta = lambd + mu

delta = 0.01

xvals = np.linspace(0, Lx, num=Nx)
yvals = np.linspace(0, Ly, num=Ny)
zvals = np.linspace(0, Lz, num=Nz)

surface = np.zeros((Nx, Ny))
sol = np.zeros((Nx, Ny, Nz))

ys, xs = np.meshgrid(xvals, yvals)
zv, yv, xv = np.meshgrid(xvals, yvals, zvals)

pressure = np.sin(xvals) + np.cos(yvals)

def ufunc(x, y, z, delta):
    parens = sigma / eta + np.square(z) / (np.square(z) + np.square(x) + np.square(y) + (delta ** 2))
    return - 1 / (4 * math.pi * mu * (np.square(x) + np.square(y) + np.square(z) + (delta ** 2))) * parens

def vfunc(x, y, z, delta):
    parens = 1 + z / (np.square(x) + np.square(y) + np.square(z) + (delta ** 2)) + eta * (np.square(x) + np.square(y)) * z / (mu * np.power(np.square(x) + np.square(y) + np.square(z) + (delta ** 2), 3))
    return - 1 / (4 * math.pi * eta * np.sqrt(np.square(x) + np.square(y) + (delta ** 2))) * parens

def sfunc(x, y, z, delta):
    xpart = x * vfunc(x, y, z, delta) / np.sqrt(np.square(x) + np.square(y) + (delta ** 2))
    ypart = y * vfunc(x, y, z, delta) / np.sqrt(np.square(x) + np.square(y) + (delta ** 2))
    return xpart, ypart, ufunc(x, y, z, delta)

def Gfunc(x, xd, y, yd, z, delta):
    return sfunc(x - xd, y - yd, z, delta)

for i in range(Nx):
    for j in range(Ny):
        for k in range(Nz):
            print(i, j, k)
            sol[i, j, k] = np.sum(Gfunc(xvals[i], xv, yvals[j], yv, zvals[k], delta) * np.moveaxis(np.tile(pressure, (Nz, 1, 1)), 0, 2) * dx * dy)

plt.contourf(sol[i, j, 0])
