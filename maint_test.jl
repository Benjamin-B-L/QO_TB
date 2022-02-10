using LinearAlgebra
using BenchmarkTools

include("params.jl")
include("layers.jl")
include("bareFS.jl")
include("dos_wmagfield.jl")

## Simulation Parameters
nx         = 100                  #System size in x (and y) direction
bmin       = 0                    #Minimum magnetic field (T)
bmax       = 100                  #Maximum magnetic field (T)
nb         = 100                  #Nbr of field points
theta_in   = 60                   #Field incident angle in degree
eta        = 0.005                 #Broadening
outfolder  = "test"               #Folder where output is written
nk         = 100                     #Nbr of kpoints for the bare FS calculation.

## Layers parameters
nlayer     = 2                    #Nbr of layers
t0         = [1.0,1.0]                #Bare NN hopping
tp0        = [-0.2,-0.2]                #Bare NNN hopping
J          = [0.1,0.1]                #Exchange parameter (for hopping renormalization)
chi        = [0.338,0.338]                #Susceptibility (for hopping renormalization)
x          = [0.1,0.08]                #Doping (for hopping renormalization)
D          = [0.03,0.12]                #Coupling to auxiliary fermions (YRZ phenomenological model)
yrz_mode   = ["Norm","AFM"]       #Mode for YRZ model
mu         = [-0.14,-0.08]                #Chemical potential of physical electrons
mu_aux     = [-0.08,-0.08]                #Chemical potential of auxiliary electrons
Q          = Dict(1=>[[2*pi/4,0.0],[0.0,2*pi/4]],2=>[])          #CDW wave vector, for the moment implemented as bidirectional only > provide [Qx,Qy]
Pcdw       = [0.03,0.0]            #CDW coupling constant
t_orth     = [0.0]                   #Interlayer hopping, empty if nlayer=1
ac         = []                   #Interlayer distance, empty if nlayer=1


## Generate model & simulation
pSimulation = SimulationParameters(nx=nx,bmin=bmin,bmax=bmax,nb=nb,theta_in=theta_in,eta=eta,outfolder=outfolder)
pLayer = LayerParameters(nlayer=nlayer,t0=t0,tp0=tp0,Q=Q,Pcdw=Pcdw,D=D,yrz_mode=yrz_mode,mu=mu,mu_aux=mu_aux,x=x,chi=chi,J=J)

#Compute the bare Fermi Surface & write it to file
# FS = getBareFS(pLayer,nk,pSimulation.eta)
#
# f = open("test.dat","w")
# for i=1:nk
#     for j=1:nk
#         write(f,string(FS.kgrid[i,j,1])*" ")
#         write(f,string(FS.kgrid[i,j,2])*" ")
#         for orb=1:FS.norb
#             write(f,string(FS.fs[i,j,orb])*" ")
#         end
#         write(f,"\n")
#     end
# end

#Compute DOS for all magnetic field values
dos = getBdepDOS(pSimulation,pLayer)
