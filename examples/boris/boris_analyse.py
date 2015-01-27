import math
import csv
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


colors = ['b', 'r', 'y', 'c', 'm']


def distance(a, b):
    d = a - b
    return math.sqrt(np.power(d, 2).sum())


class BorisData(object):
    def __init__(self, build_dir, nsteps, niter, nparticles):
        self._build_dir = build_dir
        self._nsteps = nsteps
        self._niter = niter
        self._nparticles = nparticles

        self._x = np.ndarray((self._nparticles, self._nsteps + 1), dtype=np.double)
        self._y = np.ndarray((self._nparticles, self._nsteps + 1), dtype=np.double)
        self._z = np.ndarray((self._nparticles, self._nsteps + 1), dtype=np.double)
        self._x_c = np.ndarray(self._nsteps + 1, dtype=np.double)
        self._y_c = np.ndarray(self._nsteps + 1, dtype=np.double)
        self._z_c = np.ndarray(self._nsteps + 1, dtype=np.double)
        self._x_ref = np.ndarray(self._nsteps + 1, dtype=np.double)
        self._y_ref = np.ndarray(self._nsteps + 1, dtype=np.double)
        self._z_ref = np.ndarray(self._nsteps + 1, dtype=np.double)
        self._residual = np.ndarray(self._nsteps + 1, dtype=np.double)
        self._energy = np.ndarray(self._nsteps + 1, dtype=np.double)
        self._drift = np.ndarray(self._nsteps + 1, dtype=np.double)
        self._distances = np.ndarray((self._nparticles, self._nsteps + 1), dtype=np.double)
        self._dist_min = np.ndarray(self._nsteps + 1, dtype=np.double)
        self._dist_max = np.ndarray(self._nsteps + 1, dtype=np.double)
        self._dist_mean = np.ndarray(self._nsteps + 1, dtype=np.double)

        self._read_data()

    @property
    def nsteps(self):
        return self._nsteps

    @property
    def niter(self):
        return self._niter

    @property
    def nparticles(self):
        return self._nparticles

    @property
    def x(self):
        return self._x

    @property
    def y(self):
        return self._y

    @property
    def z(self):
        return self._z

    @property
    def x_c(self):
        return self._x_c

    @property
    def y_c(self):
        return self._y_c

    @property
    def z_c(self):
        return self._z_c

    @property
    def x_ref(self):
        return self._x_ref

    @property
    def y_ref(self):
        return self._y_ref

    @property
    def z_ref(self):
        return self._z_ref

    @property
    def residual(self):
        return self._residual

    @property
    def energy(self):
        return self._energy

    @property
    def drift(self):
        return self._drift

    @property
    def distances(self):
        return self._distances

    @property
    def dist_min(self):
        return self._dist_min

    @property
    def dist_max(self):
        return self._dist_max

    @property
    def dist_mean(self):
        return self._dist_mean

    def _get_coords(self, step, p):
        return np.asarray([self.x[p][step], self.y[p][step], self.z[p][step]])

    def _read_data(self):
        data_file = self._build_dir + "/s%d_i%d_dt0.015625_m5_p%d.csv" % (self.nsteps, self.niter, self.nparticles)
        ref_file = self._build_dir + "/s%d_i%d_dt0.015625_m5_p1.csv" % (self.nsteps, self.niter)

        with open(data_file) as csvfile:
            csv_file = csv.reader(csvfile, delimiter=',')
            for row in csv_file:
                step = int(row[0])
                iter = int(row[1])
                part = int(row[2])
                if iter == 0 and step == 1:
                    self._x[part][0] = row[3]
                    self._y[part][0] = row[4]
                    self._z[part][0] = row[5]
                    self._energy[0] = row[9]
                    self._drift[0] = row[10]
                    self._residual[0] = row[11]
                    if part == -1:
                        self._x_c[0] = row[3]
                        self._y_c[0] = row[4]
                        self._z_c[0] = row[5]
                if iter == self.niter - 1:
                    self._energy[step] = row[9]
                    self._drift[step] = row[10]
                    self._residual[step] = row[11]
                    if part == -1:
                        self._x_c[step] = row[3]
                        self._y_c[step] = row[4]
                        self._z_c[step] = row[5]
                    else:
                        self._x[part][step] = row[3]
                        self._y[part][step] = row[4]
                        self._z[part][step] = row[5]

        with open(ref_file) as csvfile:
            csv_file = csv.reader(csvfile, delimiter=',')
            for row in csv_file:
                step = int(row[0])
                iter = int(row[1])
                if iter == self.niter - 1:
                    self._x_ref[step] = row[3]
                    self._y_ref[step] = row[4]
                    self._z_ref[step] = row[5]
                if iter == 0 and step == 1:
                    self._x_ref[0] = row[3]
                    self._y_ref[0] = row[4]
                    self._z_ref[0] = row[5]

        for step in range(self._nsteps + 1):
            center = np.asarray([self.x_c[step], self.y_c[step], self.z_c[step]])
            for p in range(self._nparticles):
                self._distances[p][step] = distance(self._get_coords(step, p), center)
        self._dist_min = self._distances.min(axis=0)
        self._dist_max = self._distances.max(axis=0)
        self._dist_mean = self._distances.mean(axis=0)


def get_meetup_after(data, after=0):
    _dist = data.dist_max - data.dist_min
    return np.where(_dist[after:] == _dist[after:].min())[0][0] + after


def plot_trajectories(data, until=None):
    if until is None:
        until = data.nsteps
    ndisplay = until + 1

    print("Initial Center of Mass: % 10.4f\t% 10.4f\t% 10.4f" % (data.x_c[:1][0], data.y_c[:1][0], data.z_c[:1][0]))
    print("Last Center of Mass:    % 10.4f\t% 10.4f\t% 10.4f" % (data.x_c[until:ndisplay][0], data.y_c[until:ndisplay][0], data.z_c[until:ndisplay][0]))
    print("Initial Energy: %16.4f" % data.energy[0])
    print("Final Energy:   %16.4f" % data.energy[-1])
    print("Final Drift:    %16.4f" % data.drift[-1])
    print("Final Residual: %16.4e (max: %.4e)" % (data.residual[-1], data.residual.max()))

    plt.figure(figsize=(15, 15), dpi=1200)
    ax = plt.axes(projection='3d')
    plt.title("Trajectory of Particles\ndashed: single; black solid: reference; green solid: center of mass; dot: initial point; cross: end point")
    ax.plot(data.x_ref[:ndisplay], data.y_ref[:ndisplay], data.z_ref[:ndisplay], '--k')
    if data.nparticles <= 5:
        for p in range(data.nparticles):
            ax.plot(data.x[p][:ndisplay], data.y[p][:ndisplay], data.z[p][:ndisplay], '--%s' % colors[p])
    ax.plot(data.x_c[:ndisplay], data.y_c[:ndisplay], data.z_c[:ndisplay], '-g')
    ax.plot(data.x_c[:1], data.y_c[:1], data.z_c[:1], 'ok')
    ax.plot(data.x_c[until:ndisplay], data.y_c[until:ndisplay], data.z_c[until:ndisplay], 'xk', markeredgewidth=2, markersize=5)


def plot_analytics(data, start=0, until=None, wo_drift=False):
    if until is None:
        until = data.nsteps
    ndisplay = until + 1
    figid = 1
    maxfig = 4 if wo_drift else 5
    figheight = 10 if wo_drift else 15
    steprange = range(start, ndisplay)

    plt.figure(figsize=(15, figheight), dpi=1200)
    plt.subplots_adjust(wspace=0.1)

    plt.subplot(maxfig, 1, figid)
    plt.plot(steprange, data.energy[start:ndisplay])
    plt.grid()
    plt.xlim(start, until)
    plt.ylabel("Energy")
    figid += 1

    if wo_drift == False:
        plt.subplot(maxfig, 1, figid)
        plt.plot(steprange, data.drift[start:ndisplay] / data.energy[start:ndisplay])
        plt.grid()
        plt.xlim(start, until)
        plt.ylabel("Relative Drift\n(regarding total energy)")
        figid += 1

    plt.subplot(maxfig, 1, figid)
    plt.plot(steprange, data.residual[start:ndisplay])
    plt.grid()
    plt.xlim(start, until)
    plt.yscale("log")
    plt.ylabel("Residual")
    plt.xlabel("Step")
    figid += 1

    plt.subplot(maxfig, 1, figid)
    plt.plot(steprange, data.distances.transpose()[start:ndisplay])
    plt.grid()
    plt.xlim(start, until)
    plt.ylabel("Distance of Particles\nfrom Center of Mass")
    figid += 1

    plt.subplot(maxfig, 1, figid)
    plt.plot(steprange, data.dist_min.transpose()[start:ndisplay], '-r')
    plt.plot(steprange, data.dist_max.transpose()[start:ndisplay], '-b')
    plt.plot(steprange, data.dist_mean.transpose()[start:ndisplay], '-g')
    plt.grid()
    plt.xlim(start, until)
    plt.ylabel("Stats of Distance of Particles\nfrom Center of Mass")
    figid += 1


def plot_particle_meetup(data, step, width=10):
    start = step - width
    stop = step + width + 1

    plt.figure(figsize=(15, 15), dpi=1200)
    ax = plt.axes(projection='3d')
    plt.title("Particle Meetup at step %d\ndot: start; square: end" % step)
    if data.nparticles <= 5:
        for p in range(data.nparticles):
            ax.plot(data.x[p][start:stop], data.y[p][start:stop], data.z[p][start:stop], '--%s' % colors[p])
    ax.plot(data.x_c[start:stop], data.y_c[start:stop], data.z_c[start:stop], '-g')
    ax.plot(data.x_c[step:step+1], data.y_c[step:step+1], data.z_c[step:step+1], 'xk', markeredgewidth=3, markersize=10)
    ax.plot(data.x_c[start:start+1], data.y_c[start:start+1], data.z_c[start:start+1], 'ok')
    ax.plot(data.x_c[stop-1:stop], data.y_c[stop-1:stop], data.z_c[stop-1:stop], 'sk')

    plot_analytics(data, start, stop-1, True)
