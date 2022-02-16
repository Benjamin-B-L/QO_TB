using LinearAlgebra
using BenchmarkTools

include("params.jl")
include("layers.jl")
include("bareFS.jl")
include("dos_wmagfield.jl")

## Simulation Parameters
nx         = 1000                  #System size in x direction
ny_avg     = 5                      #Nbr of ky points for averaging
bmin       = 1/660                   #Minimum magnetic field (T)
bmax       = 1/580                  #Maximum magnetic field (T)
nb         = 5000                  #Nbr of field points
theta_in   = 0.0                   #Field incident angle in degree
eta        = 0.0001                 #Broadening
outfolder  = "test"               #Folder where output is written
nk         = 100                     #Nbr of kpoints for the bare FS calculation.

## Layers parameters
nlayer     = 1                    #Nbr of layers
t0         = [1.0]                #Bare NN hopping
tp0        = [0.0]                #Bare NNN hopping
J          = [0.0]                #Exchange parameter (for hopping renormalization)
chi        = [0.0]                #Susceptibility (for hopping renormalization)
x          = [1.0]                #Doping (for hopping renormalization)
D          = [0.0]                #Coupling to auxiliary fermions (YRZ phenomenological model)
yrz_mode   = ["AFM"]       #Mode for YRZ model
mu         = [-0.2]                #Chemical potential of physical electrons
mu_aux     = [-0.5]                #Chemical potential of auxiliary electrons
Q          = Dict(1=>[])          #CDW wave vector, for the moment implemented as bidirectional only > provide [Qx,Qy]
Pcdw       = [0.0]            #CDW coupling constant
g_tilde    = [0.0]            #Zeeman splitting constant
t_orth     = []                   #Interlayer hopping, empty if nlayer=1
ac         = []                   #Interlayer distance, empty if nlayer=1


## Generate model & simulation
pSimulation = SimulationParameters(nx=nx,ny_avg=ny_avg,bmin=bmin,bmax=bmax,nb=nb,theta_in=theta_in,eta=eta,outfolder=outfolder)
pLayer = LayerParameters(nlayer=nlayer,t0=t0,tp0=tp0,Q=Q,Pcdw=Pcdw,D=D,yrz_mode=yrz_mode,mu=mu,mu_aux=mu_aux,x=x,chi=chi,J=J,t_orth=t_orth,ac=ac)



# #Compute the bare Fermi Surface & write it to file
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

# #Compute DOS for all magnetic field values
# dos = getBdepDOS(pSimulation,pLayer)
#
#
# f = open("test_dos_nx1000_nb5000_maharajparams.dat","w")
# for ib = 1:nb
#     write(f,string(1/pSimulation.bgrid[ib])*" ")
#     for iy = 1:ny_avg+1
#         write(f,string(dos.dos[ib,iy])*" ")
#     end
#     write(f,"\n")
# end
# close(f)
