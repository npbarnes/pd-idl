import numpy as np
from scipy.integrate import ode,quad

class mom_trans:

    def __init__(self,  rp=1100, nsw=0.0273e15, mp=1.67e-27, msw=None, mpu=None, vsw=400, 
                        B0=0.3e-9, Ly=None, max_ndot=5e11,cap_r=None):
        self.rp = rp
        self.nsw = nsw
        self.mp = mp
        self.msw = msw if msw is not None else 1*self.mp               # Mass of average solar wind particle
        self.mpu = mpu if mpu is not None else 16*self.mp              # Mass of average pick up ion
        self.vsw = vsw                                                 # Solar wind speed
        self.B0 = B0                                                   # Ambient magnetic field in T
        self.Ly = Ly if Ly is not None else 2*self.rp                  # Field aligned extend of the cloud
        self.max_ndot = max_ndot                                       # Peak number of ions per second picked up
        self.cap_r = 1.05*self.rp

        self.rho = self.msw*self.nsw
        self.vA = (self.B0/np.sqrt(np.pi*4e-7*self.rho/1e9))/1e3 # Alven speed in km/s

    def rho_pu_dot(self, r):
        return self.mpu*self.max_ndot*np.exp(-r**2/(3*self.rp)**2)

    def rho_pu_dot2(self, r):
        r = max(r, self.cap_r)
        return 1e-15*self.mpu*self.max_ndot*(1e15*(self.rp/r)**25 + 5e9*(self.rp/r)**8)
        

    def f(self, t, s):
        x = s[0]
        v = s[1]
        rho_pu = s[2]

        r = np.sqrt(x**2)

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


    def _plot_energies(self, ax, s_arr, mass_list, kwargs_list):
        for mass, kwargs in zip(mass_list, kwargs_list):
            ax.plot(-s_arr[:,0]/self.rp, 0.5*mass*(s_arr[:,1]*1e3)**2/1.6e-19, **kwargs)

    def plot_energies(self, ax, mass_list, kwargs_list):
        t_arr, s_arr = self.integrate_until(10000,200*self.rp)
        self._plot_energies(ax, s_arr, mass_list, kwargs_list)

