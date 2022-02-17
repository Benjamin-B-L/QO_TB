# if core_QOparallel is in another folder
# import sys
# sys.path.insert(1,'$path')
from core_QOparallel import SimulationParameters, LayerParameters, call_parallelQO

# TYPE OF CALCULATION
BareFSswitch = False
DOSswitch = True
# SIMULATION PARAMETERS
nx = 1000                           # Number of site in x direction
ny_avg = 5                          # Number of ky points for averaging
bmin = 1/660                        # Minimum magnetic field value
bmax = 1/580                        # Maximum magnetic field value
nb = 5000                           # Number of magnetic field points
theta_in = 0.0                      # Angle of the magnetic field (degree)
eta = 0.0001                        # Broadening
nk = 100                            # Nbr of kpoints for bare FS
outfile = 'test'                    # Output file

# LAYER PARAMETERS
nlayer = 1                          # Number of layers
t0 = [1.0]                          # NN hopping
tp0 = [0.0]                         # NNN hopping
J = [0.0]                           # Spin-exchange for hopping renorm
chi = [0.0]                         # Susceptibility for hopping renorm
x = [1.0]                           # Doping for hopping renorm
D = [0.0]                           # YRZ coupling to auxilary fermions
yrz_mode = ["NORM"]                 # NORM or AFM type of YRZ coupling
mu = [-0.2]                         # Chemical potential physical fermions
mu_aux = [0.0]                      # Chemical potential auxiliary fermions
Q = {1: []}                         # CDW vector
Pcdw = [0.0]                        # CDW field
g_tilde = [0.0]                     # Zeeman splitting factor
t_orth = []                         # interlayer coupling
ac = []                             # interlayer distance c/a
# Number of cores
ncores = 4


# Construct parameters structs
pSim = SimulationParameters(nx=nx, ny_avg=ny_avg, bmin=bmin, bmax=bmax,
                            nb=nb, theta_in=theta_in, eta=eta, nk=nk, outfile=outfile)
pLayer = LayerParameters(nlayer=nlayer, t0=t0, tp0=tp0, J=J, chi=chi, x=x, D=D,
                         yrz_mode=yrz_mode, mu=mu, mu_aux=mu_aux, Q=Q, Pcdw=Pcdw,
                         g_tilde=g_tilde, t_orth=t_orth, ac=ac)
# Compute DOS in parallel (bgrid is sliced and distributed over ncores)
if __name__ == "__main__":
    call_parallelQO(DOSswitch, BareFSswitch, pSim, pLayer, ncores)
