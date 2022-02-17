###Type defs
const FieldMesh = Array{Float64,1}

########## Parameters struct ###########

"""
    SimulationParameters
    Contains the general parameters necessary for the simulation

Fields:
----------

- nx              : System size in x-direction (sets also the nbr of k-points for FS calculation and ky random determination)
- ny              : Nbr of ky-points for averaging DOS
- bmin            : Minimun value for the external magnetic field
- bmax            : Maximum value for the external magnetic field
- nb              : Number of points between bmin & bmax
- theta           : Magnetic field polar angle (give in Degree, is then converted to rad)
- eta             : broadening parameter (inverse of particle lifetime)
- nk              : Number of k-points for bare FS
- outfolder       : Folder where output files will be written
- bgrid           : Magnetic field array [bmin,bmax] (Float64,nb)
"""

struct SimulationParameters
    nx          :: Int64
    ny_avg      :: Int64
    bmin        :: Float64
    bmax        :: Float64
    nb          :: Int64
    theta       :: Float64
    eta         :: Float64
    nk          :: Int64
    outfile   :: AbstractString
    bgrid       :: FieldMesh

    SimulationParameters(;nx,ny_avg,bmin,bmax,nb,theta_in,eta,nk,outfile) = nx > 0 && ny_avg > 0 && nb > 0 ? new(nx,ny_avg,bmin,bmax,
                                                                                            nb,theta_in*pi/180,
                                                                                            eta,nk,outfile,
                                                                                            [bmin+n*(bmax-bmin)/nb for n=0:nb-1]) :
                        throw(ArgumentError("nx=$nx (system size), ny_avg=$ny_avg (nbr of averaging points) and nb=$nb (Nbr of field points) should be larger than zero !"))
end
