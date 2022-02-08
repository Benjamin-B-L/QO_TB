include("params.jl")
include("layers.jl")

## Simulation Parameters
nx         = 100                  #System size in x (and y) direction
bmin       = 0                    #Minimum magnetic field (T)
bmax       = 100                  #Maximum magnetic field (T)
nb         = 1000                 #Nbr of field points
theta_in   = 60                   #Field incident angle in degree
eta        = 0.01                 #Broadening
outfolder  = "test"               #Folder where output is written
Nk         = 100                  #Nbr of kpoints for the bare FS calculation.

## Layers parameters
nlayer     = 1                    #Nbr of layers
t0         = [1.0]                #Bare NN hopping
tp0        = [0.0]                #Bare NNN hopping
J          = [0.0]                #Exchange parameter (for hopping renormalization)
chi        = [0.0]                #Susceptibility (for hopping renormalization)
x          = [1.0]                #Doping (for hopping renormalization)
D          = [0.0]                #Coupling to auxiliary fermions (YRZ phenomenological model)
mu         = [0.0]                #Chemical potential of physical electrons
mu_aux     = [0.0]                #Chemical potential of auxiliary electrons
Q          = Dict(1=>[])          #CDW wave vector, for the moment implemented as bidirectional only > provide [Qx,Qy]
t_orth     = []                   #Interlayer hopping, empty if nlayer=1
ac         = []                   #Interlayer distance, empty if nlayer=1


## Generate model & simulation
pSimulation = SimulationParameters(nx=nx,bmin=bmin,bmax=bmax,nb=nb,theta_in=theta_in,eta=eta,outfolder=outfolder)
pLayer = LayerParameters(nlayer=nlayer,t0=t0,tp0=tp0)
