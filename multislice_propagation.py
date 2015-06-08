# I copied this from for_crif.py
from time import *
from pyspecdata import *
from pyspecdata.prop import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
I_z = 0.5*array([[1,0],[0,-1]],dtype='complex128')
E = array([[1,0],[0,1]],dtype='complex128')
I_x = 0.5*array([[0,1],[1,0]],dtype='complex128')
I_y = 0.5*array([[0,1j],[-1j,0]],dtype='complex128')
I_p = array([[0,1],[0,0]],dtype='complex128')
# note that i can test fortran order with .flags, and that it has no effect on 1D arrays!
fl = figlist()

n_nu = 200

# do an adiabatic pulse
beta = 5.4
multiplier = 1/1.5167e-3*100/234*17
scaling_factor = 1.
scaleup_freq = 2
nu_sweep = 100e6
nu_plot = 1.5*nu_sweep
B1_inv = 48e6/scaling_factor
B1_exc = 9.6/48.*B1_inv
tau_pulse = 620e-9*(scaling_factor**2)
#(28.1175e9/1e4)
print 'nu_sweep',nu_sweep, 'B1_exc',B1_exc, 'B1_inv',B1_inv, 'tau_pulse',tau_pulse
t_steps = 500
nu_offset = linspace(-nu_plot,nu_plot,n_nu)
#t = linspace(0,2.5*tau_pulse,2.5*t_steps)
t = r_[0.:3*tau_pulse:1e-9]
#t_axis = linspace(0,tau_pulse*3,t_steps*3)
t_axis = t.copy()
fpg = (28.1175e9/1e4)
mu = -2*pi*nu_sweep/beta
bprime = beta*2./tau_pulse
complex_form = (sech(bprime * (t-tau_pulse/2)))**(1.+1j*mu) # no idea why this isn't working
dt = t[1]-t[0]
#{{{ generate the complex for the second
pulse_center = tau_pulse/2
tau = (t-pulse_center)*2/tau_pulse
fmod = nu_sweep*tanh(beta*tau)
fmod[t<pulse_center - tau_pulse/2] = 0
fmod[t>pulse_center + tau_pulse/2] = 0
pulse_phase = dt * cumsum(fmod)# second pulse
pulse_phase = exp(1j*2.*pi*pulse_phase)
complex_form1 = B1_exc * sech(beta*tau) * pulse_phase
complex_form1[t<pulse_center - tau_pulse/2] = 0
complex_form1[t>pulse_center + tau_pulse/2] = 0
#}}}
fmod1 = fmod
#{{{ generate the complex for the second
pulse_center = 3*tau_pulse/2
tau = (t-pulse_center)*2/tau_pulse
fmod = nu_sweep*tanh(beta*tau)
fmod[t<pulse_center - tau_pulse/2] = 0
fmod[t>pulse_center + tau_pulse/2] = 0
pulse_phase = dt * cumsum(fmod)# second pulse
pulse_phase = exp(1j*2.*pi*pulse_phase)
complex_form2 = B1_inv * sech(beta*tau) * pulse_phase
complex_form2[t<pulse_center - tau_pulse/2] = 0
complex_form2[t>pulse_center + tau_pulse/2] = 0
#}}}
complex_form = complex_form1 + complex_form2 # the overlay of the two pulses
fl.next('waveforms')
plot(abs(complex_form)/fpg,'k',linewidth=3)
plot(abs(complex_form1)/fpg,'r',linewidth=1.5)
plot(abs(complex_form2)/fpg,'b',linewidth=1.5)
fl.next('coherence')
t_dep = array([r_[(fmod+fmod1)] , r_[abs(complex_form)] , r_[ones(shape(t))]]) * tau_pulse / t_steps
operators = array([I_z          , I_x                   , I_z])
parameters = (array([1.])       , array([1.])           , nu_offset) # these need to have the same number of dimensions as the full parameter grid

# generate the density matrices, w/out time propagation
#{{{ excitation
U = prop(t = t_dep,p = parameters,H = operators,accumtime = True)
rho = sandwich(U[:,:,:,:],
	array(I_z,dtype='complex128')
	)
#}}}
signal = lvdot(rho,I_p)
signal_mask = abs(signal)
image(abs(signal),
		x=t_axis/1e-9,y=nu_offset/1e6)
xlabel(r't / $ns$')
ylabel(r'$\Omega$ / 2 $\pi\;MHz$')
title('Magnitude of coherence')
fl.next('signal_complex')
image(signal,
		x=t_axis/1e-9,y=nu_offset/1e6)
xlabel(r't / $ns$')
ylabel(r'$\Omega$ / 2 $\pi\;MHz$')
title('complex coherence')
fl.next('inversion')
image(abs(lvdot(rho,I_z)),
		x=t_axis/1e-9,y=nu_offset/1e6)
xlabel(r't / $ns$')
ylabel(r'$\Omega$ / 2 $\pi\;MHz$')
title('inversion')
fl.next('signal')
plotfunc = angle(lvdot(rho,I_p))
plotfunc = diff(plotfunc,axis=0)**2
plotfunc = -log(abs(plotfunc))
plotfunc *= signal_mask[:-1,:]
print plotfunc
nu_offset = nu_offset[:-1]
print shape(plotfunc)
image(plotfunc,
		x=t_axis/1e-9,y=nu_offset/1e6)
xlabel(r't / $ns$')
ylabel(r'$\Omega$ / 2 $\pi\;MHz$')

print 'nu_sweep',nu_sweep/1e6, 'MHz, B1_exc',B1_exc/1e6, 'MHz B1_inv',B1_inv/1e6,'MHz --> ',B1_inv/(28.1175e9/1e4),'G tau_pulse',tau_pulse/1e-9,'ns dt',dt*1e9,'ns'
title('Signal Level (log of inverse of dispersion)')
fl.show('prop_test_130215')
