import numpy as np
from scipy.integrate import ode

class mom_trans:

    def __init__(self, rp=None, nsw=None, mp=None, msw=None, mpu=None, vsw=None, B0=None, Ly=None, max_ndot=None, y=None):
        self.rp = 1100. if rp is None else rp                     # Pluto radius
        self.nsw = 0.0273*1e15 if nsw is None else nsw            # number density of the solar wind in km^-3
        self.mp = 1.67e-27 if mp is None else mp                  # mass of a proton
        self.msw = 1*self.mp if msw is  None else msw                  # Mass of average solar wind particle
        self.mpu = 16*self.mp if mpu is None else mpu                  # Mass of average pick up ion
        self.vsw = 400.0 if vsw is None else vsw                  # Solar wind speed
        self.B0 = 0.3e-9 if B0 is None else B0                    # Ambient magnetic field in T
        self.Ly = 2.0*self.rp if Ly is None else Ly                    # Field aligned extend of the cloud
        self.max_ndot = 5e11 if max_ndot is None else max_ndot    # Peak number of ions per second picked up
        self.y = 0 if y is None else y                            # Offset from the centerline

        self.rho = self.msw*self.nsw
        self.vA = (self.B0/np.sqrt(np.pi*4e-7*self.rho/1e9))/1e3 # Alven speed in km/s

    def rho_pu_dot(self, r):
        return self.mpu*self.max_ndot*np.exp(-r**2/(3*self.rp)**2)
        

    def f(self, t, s):
        x = s[0]
        v = s[1]
        rho_pu = s[2]

        r = np.sqrt(x**2 + self.y**2)

        return np.array([v,
                         ((self.vsw - v)*self.vA*self.rho/((self.rho+rho_pu)*self.Ly) 
                             - self.rho_pu_dot(r)*v/(self.rho+rho_pu)),
                         self.rho_pu_dot(r)])

    def integrate_until(self, t, x):
        r = ode(self.f).set_integrator('dopri5')
        r.set_initial_value(np.array([-30*self.rp, self.vsw, 0]), 0)

        t_list = []
        s_list = []
        while(r.successful() and r.t < t and r.y[0]<x):
            t_list.append(r.t)
            s_list.append(r.y)
            r.integrate(r.t+1)

        t_arr = np.array(t_list)
        s_arr = np.array(s_list)

        return t_arr, s_arr

    def traj_energies(self, ax, x1, y1, x2, y2, mass_list, kwargs_list):
        slp = (y2-y1)/(x2-x1)
        xlst = np.linspace(-20*self.rp, 200*self.rp, 500)
        ylst = -slp*xlst + y1*self.rp + 25*self.rp

        slst = []
        for x, y in zip(xlst, ylst):
            self.y = y
            t, s = self.integrate_until(10000, x)
            slst.append(s[-1,:])

        slst = np.array(slst)

        self._plot_energies(ax, slst, mass_list, kwargs_list)


    def _plot_energies(self, ax, s_arr, mass_list, kwargs_list):
        for mass, kwargs in zip(mass_list, kwargs_list):
            ax.plot(-s_arr[:,0]/self.rp, 0.5*mass*(s_arr[:,1]*1e3)**2/1.6e-19, **kwargs)

    def plot_energies(self, ax, mass_list, kwargs_list):
        t_arr, s_arr = self.integrate_until(10000,200*self.rp)
        self._plot_energies(ax, s_arr, mass_list, kwargs_list)

