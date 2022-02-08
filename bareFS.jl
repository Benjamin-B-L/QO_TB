include("layers.jl")

###Type defs
const FSArray = Array{Float64,3}
const OrbsArray = Array{AbstractString,1}

########## FS struct ###########

"""
Fermi Surface
Contains Fermi surface data calculated without magnetic field

Fields:
----------

- nk              : Nbr of k-points per direction
- norb            : Nbr of orbitals
- kgrid           : k-point grid (set up for 2D square BZ)
- fs              : Fermi surface

"""

struct FermiSurface

    nk        :: Int64
    norb      :: Int64
    ncdw      :: Int64
    kgrid     :: FSArray
    orbs      :: OrbsArray
    fs        :: FSArray

    function FermiSurface(;nk,norb,ncdw,orbs)
        # Check nbr of k-points
        if nk <= 0 || norb <= 0
            throw(ArgumentError("nk=$nk number of kpoints and norb=$norb should be larger than zero"))
        end
        kgrid = zeros(Float64,nk,nk,2)
        for ikx = 1:nk
            for iky = 1:nk
                kgrid[ikx,iky,1] = -pi+2*pi*(ikx-1)/nk
                kgrid[ikx,iky,2] = -pi+2*pi*(iky-1)/nk
            end
        end
        return new(nk,norb,ncdw,kgrid,orbs,zeros(Float64,nk,nk,norb))
    end
end


################# Functions #################################################

"""
    getBareFS(Layer::LayerParam,nk::Int64)

Calculates the Fermi Surface without magnetic field from the Layers' parameters
Layer, the number of k-points nk.
"""
function getBareFS(Layer::LayerParameters,nk::Int64,eta::Float64)
    #Initialize FS struct
    norb,ncdw,orbs = prepFSHam(Layer)
    FS = FermiSurface(nk=nk,norb=norb,ncdw=ncdw,orbs=orbs)

    #Compute FS at each k-point
    for ikx = 1:nk
        for iky = 1:nk
            FS.fs[ikx,iky,:] = getRhoatmu(Layer,norb,ncdw,FS.kgrid[ikx,iky,:],eta)
        end
    end
    return FS
end

"""
    prepFSHam(Layer::LayerParam)

Determines the size of the Hamiltonian in k-space for FS calculation
"""
function prepFSHam(Layer::LayerParameters)
    #First search for CDW
    tmp = []
    for ilayer = 1:nlayer
        if length(Layer.Q[ilayer]) > 0
            for iq = 1:length(Layer.Q[ilayer])
                append!(tmp,Int64(2*pi/sum(Layer.Q[ilayer][iq])))
            end
        end
    end
    if length(tmp)>0
        ncdw = prod(tmp)
    else
        ncdw = 0
    end

    #Compute norb by checking if YRZ & construct orbs array
    norb=0
    orbs = String[]
    orb_cnt=1
    for ilayer = 1:nlayer
        if Layer.D[ilayer] < 1e-8
            norb += max(ncdw,1)              #No auxiliary fermions in layer ilayer
            for icdw = 1:max(ncdw,1)         #get all orbs
                append!(orbs,["c"*string(ilayer)*"_k"*string(icdw)])
            end
        else
            norb += 2*max(ncdw,1)            #Auxiliary fermions in layer ilayer
            for icdw = 1:max(ncdw,1)
                append!(orbs,["c"*string(ilayer)*"_k"*string(icdw)])
            end
            for icdw = 1:max(ncdw,1)
                append!(orbs,["ctilde"*string(ilayer)*"_k"*string(icdw)])
            end
        end
    end

    return norb,ncdw,orbs
end

"""
    getRhoatmu(Layer::LayerParameters,norb,ncdw,k)

Compute the spectral weight at the fermi level at k point from the layers' parameters in Layer. The
size of Hamiltonian is given by norb, and there is (2*)ncdw k-points per layer.
"""
function getRhoatmu(Layer::LayerParameters,norb::Int64,ncdw::Int64,k::Vector{Float64}, eta::Float64)
    kx, ky = k
    #evaluate Green's function at mu
    rho = inv(1im*eta*I(norb))
    println(rho)
    #
    return zeros(Float64,norb)
end
