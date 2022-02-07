###Type defs
const FieldMesh = Array{Float64,1}

########## Parameters struct ###########

"""
    SimulationParameters
    Contains the general parameters necessary for the simulation

Fields:
----------

- nx              : System size in x-direction (sets also the nbr of k-points for FS calculation and ky random determination)
- bmin            : Minimun value for the external magnetic field
- bmax            : Maximum value for the external magnetic field
- nb              : Number of points between bmin & bmax
- theta           : Magnetic field polar angle (give in Degree, is then converted to rad)
- eta             : broadening parameter (inverse of particle lifetime)
- outfolder       : Folder where output files will be written
"""

struct SimulationParameters
    nx          :: Int64
    bmin        :: Float64
    bmax        :: Float64
    nb          :: Int64
    theta       :: Float64
    eta         :: Float64
    outfolder   :: AbstractString
    b           :: FieldMesh

    SimulationParameters(;nx,bmin,bmax,nb,theta_in,eta,outfolder) = nx > 0 && nb > 0 ? new(nx,bmin,bmax,
                                                                                            nb,theta_in*pi/180,
                                                                                            eta,outfolder,
                                                                                            [bmin+n*(bmax-bmin)/nb for n=0:nb-1]):
                        throw(ArgumentError("nx=$nx (system size) and nb=$nb (Nbr of field points) should be larger than zero !"))
end
