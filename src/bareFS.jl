using LinearAlgebra

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
    kgrid     :: FSArray
    orbs      :: OrbsArray
    fs        :: FSArray

    function FermiSurface(;nk,norb,orbs)
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
        return new(nk,norb,kgrid,orbs,zeros(Float64,nk,nk,norb))
    end
end


################# Functions #################################################

"""
    prepFSHam(Layer::LayerParam)

Determines the size of the Hamiltonian in k-space for FS calculation
"""
function prepFSHam(Layer::LayerParameters)
    #First search for CDW
    qlist,ncdw_list = get_qlist(Layer)

    #Get the k-vector list
    klist = get_klist(ncdw_list,qlist)
    ncdw = length(klist)

    #Compute norb by checking if YRZ & construct orbs array
    norb=0
    orbs = String[]
    orbs_layer = Int64[]
    orb_cnt=1
    for ilayer = 1:Layer.nlayer
        if Layer.D[ilayer] < 1e-8
            norb += max(ncdw,1)              #No auxiliary fermions in layer ilayer
            for icdw = 1:max(ncdw,1)         #get all orbs
                append!(orbs,["c"*string(ilayer)*"_k"*string(icdw)])
                append!(orbs_layer,ilayer)
            end
        else
            norb += 2*max(ncdw,1)            #Auxiliary fermions in layer ilayer
            for icdw = 1:max(ncdw,1)
                append!(orbs,["c"*string(ilayer)*"_k"*string(icdw)])
                append!(orbs_layer,ilayer)
            end
            for icdw = 1:max(ncdw,1)
                append!(orbs,["ctilde"*string(ilayer)*"_k"*string(icdw)])
                append!(orbs_layer,ilayer)
            end
        end
    end

    return norb,orbs,orbs_layer,klist
end

"""
    get_klist(ncdw_list::Vector{Any},qlist::Vector{Any})

Determines the list of k-points necessary for FS calculation from the q vectors
"""
function get_klist(ncdw_list::Vector{Any},qlist::Vector{Any})
    klist=[[0.0,0.0]]
    klist_round=[[0.0,0.0]]
    #If no cdw vector
    if length(qlist) == 0
        return klist
    end
    #Else, find all possible k-points
    ktmp = [[0.0,0.0]]
    iq=1
    while iq <= length(qlist)
        for k in ktmp
            for inq=1:ncdw_list[iq]
                if !(round.([mod((k+qlist[iq])[1],2pi),mod((k+qlist[iq])[2],2pi)],digits=4) in klist_round)
                    append!(klist,[[mod((k+qlist[iq])[1],2pi),mod((k+qlist[iq])[2],2pi)]])
                    append!(klist_round,[round.([mod((k+qlist[iq])[1],2pi),mod((k+qlist[iq])[2],2pi)],digits=4)])
                end
            end
        end
        #println(klist)
        if length(ktmp)==length(klist)
            if iq==length(qlist)
                break
            else
                ktmp=[klist[ik] for ik=1:length(klist)]
                iq+=1
            end
        else
            ktmp=[klist[ik] for ik=1:length(klist)]
        end
    end
    return klist
end

"""
    get_qlist(Layer::LayerParameters)

Find all the different q vectors in all layers
"""
function get_qlist(Layer::LayerParameters)
    #Extract all different Q vectors
    qlist = []
    for ilayer = 1:Layer.nlayer
        if length(Layer.Q[ilayer]) > 0
            for iq = 1:length(Layer.Q[ilayer])
                if !(Layer.Q[ilayer][iq] in qlist)
                    append!(qlist,[Layer.Q[ilayer][iq]])
                end
            end
        end
    end

    #If no cdw vector
    if length(qlist)==0
        return [],[]
    end
    #For all different Q, find the necessary nbr of kpoints
    ncdw_list = []
    for q in qlist
        if abs(q[1])>1e-8 && abs(q[2])<1e-8               #if x unidirectional
            if abs(trunc(Int64,2*pi/abs(q[1]))-2*pi/abs(q[1]))<1e-4
                append!(ncdw_list,trunc(Int64,2*pi/abs(q[1])))
            else
                throw(ArgumentError("Q vector provided is not an integer..."))
            end
        elseif abs(q[2])>1e-8 && abs(q[1])<1e-8           #if x unidirectional
            if abs(trunc(Int64,2*pi/abs(q[2]))-2*pi/abs(q[2]))<1e-4
                append!(ncdw_list,trunc(Int64,2*pi/abs(q[2])))
            else
                throw(ArgumentError("Q vector provided is not an integer..."))
            end
        elseif abs(q[1])>1e-8 && abs(q[2])>1e-8           #if (x,y) bi-directional
            if abs(q[1])-abs(q[2])>1e-8
                throw(ArgumentError("Q=(qx,qy) with abs(qx)!=abs(qy) is not yet implemented !"))
            else
                append!(ncdw_list,Int64(2*pi/abs(q[1])))
            end
        end
    end
    #Total nbr is the product
    return qlist,ncdw_list
end




"""
    eps(k,t,tp,mu)

Computes the relation dispersion at k, from NN hopping t, NNN tp and chemical potential mu
"""
function eps(k,t,tp,mu)
    return ComplexF64(-2*t*(cos(k[1])+cos(k[2]))-2*tp*(cos(k[1]+k[2])+cos(k[1]-k[2]))-mu)
end

"""
    eps_aux(k,t,tp,mu,yrz_mode)

Computes the relation dispersion at k, from NN hopping t, NNN tp and chemical potential mu
"""
function eps_aux(k,t,tp,mu,yrz_mode)
    if yrz_mode=="Norm"
        return -ComplexF64(-2*t*(cos(k[1])+cos(k[2]))-mu)
    elseif yrz_mode=="AFM"
        k+=[pi,pi]
        return ComplexF64(-2*t*(cos(k[1])+cos(k[2]))-2*tp*(cos(k[1]+k[2])+cos(k[1]-k[2]))-mu)
    else
        throw(ArgumentError("YRZ mode $yrz_mode not implemented"))
    end
end

"""
    yrz_delta(k,D)

Computes the YRZ coupling constant at k
"""
function yrz_delta(k,D)
    return D*(cos(k[1])-cos(k[2]))
end

"""
    ffactor(k)

Computes the CDW form factor (d-wave here)
"""
function ffactor(k)
    return (cos(k[1])-cos(k[2]))
end

"""
    getH(Layer::LayerParameters,norb::Int64,klist::Vector{Vector{Float64}})

Generates the Hamiltonian at zero field
"""
function getH(Layer::LayerParameters,norb::Int64,klist::Vector{Vector{Float64}})
    nklist =length(klist)
    H = zeros(Complex{Float64},norb,norb)
    cnt=1
    cnt_perlay=[]
    for ilay = 1:Layer.nlayer
        push!(cnt_perlay,cnt-1)
        #if YRZ model
        if abs(Layer.D[ilay])>1e-8
            for iyrz = 1:2
                for ik = 1:nklist
                    #Physical electrons
                    if iyrz==1
                        #Diagonal
                        H[cnt,cnt] = eps(klist[ik],Layer.t[ilay],Layer.tp[ilay],Layer.mu[ilay])
                        #YRZ coupling
                        H[cnt,cnt+nklist] = yrz_delta(klist[ik],Layer.D[ilay])
                        #CDW
                        for q in Layer.Q[ilay]
                            #+Q
                            kq = [mod((klist[ik]+q)[1]+pi,2pi)-pi,mod((klist[ik]+q)[2]+pi,2pi)-pi]
                            kq_index = cnt_perlay[ilay] + findfirst(x->abs((x[1]-kq[1])^2+(x[2]-kq[2])^2)<1e-8,klist)
                            H[cnt,kq_index] += Layer.Pcdw[ilay]*ffactor(klist[ik])
                            H[kq_index,cnt] += Layer.Pcdw[ilay]*ffactor(klist[ik])
                            #-Q
                            kq = [mod((klist[ik]-q)[1]+pi,2pi)-pi,mod((klist[ik]-q)[2]+pi,2pi)-pi]
                            kq_index = cnt_perlay[ilay] + findfirst(x->abs((x[1]-kq[1])^2+(x[2]-kq[2])^2)<1e-8,klist)
                            H[cnt,kq_index] += Layer.Pcdw[ilay]*ffactor(klist[ik])
                            H[kq_index,cnt] += Layer.Pcdw[ilay]*ffactor(klist[ik])
                        end
                        #interlayer coupling
                        if Layer.nlayer > 1
                            if ilay==1
                                H[cnt,cnt+2*nklist] = -Layer.t_orth[ilay]
                            elseif ilay==Layer.nlayer
                                H[cnt,cnt-(cnt_perlay[ilay]-cnt_perlay[ilay-1])] = -Layer.t_orth[ilay-1]
                            else
                                H[cnt,cnt+2*nklist] = -Layer.t_orth[ilay]
                                H[cnt,cnt-(cnt_perlay[ilay]-cnt_perlay[ilay-1])] = -Layer.t_orth[ilay-1]
                            end
                        end
                    #Auxiliary electrons
                    else
                        #Diagonal
                        H[cnt,cnt] = eps_aux(klist[ik],Layer.t[ilay],Layer.tp[ilay],
                                             Layer.mu_aux[ilay],Layer.yrz_mode[ilay])
                        #YRZ
                        H[cnt,cnt-nklist] = yrz_delta(klist[ik],Layer.D[ilay])
                    end
                    cnt += 1
                end
            end
        #if no YRZ
        else
            #Only physical electrons
            for ik = 1:nklist
                #Diagonal
                H[cnt,cnt] = eps(klist[ik],Layer.t[ilay],Layer.tp[ilay],Layer.mu[ilay])
                #CDW
                for q in Layer.Q[ilay]
                    #+Q
                    kq = [mod((klist[ik]+q)[1]+pi,2pi)-pi,mod((klist[ik]+q)[2]+pi,2pi)-pi]
                    kq_index = cnt_perlay[ilay] + findfirst(x->abs((x[1]-kq[1])^2+(x[2]-kq[2])^2)<1e-8,klist)
                    H[cnt,kq_index] += Layer.Pcdw[ilay]*ffactor(klist[ik])
                    H[kq_index,cnt] += Layer.Pcdw[ilay]*ffactor(klist[ik])
                    #-Q
                    kq = [mod((klist[ik]-q)[1]+pi,2pi)-pi,mod((klist[ik]-q)[2]+pi,2pi)-pi]
                    kq_index = cnt_perlay[ilay] + findfirst(x->abs((x[1]-kq[1])^2+(x[2]-kq[2])^2)<1e-8,klist)
                    H[cnt,kq_index] += Layer.Pcdw[ilay]*ffactor(klist[ik])
                    H[kq_index,cnt] += Layer.Pcdw[ilay]*ffactor(klist[ik])
                end
                #interlayer coupling
                if Layer.nlayer > 1
                    if ilay==1
                        H[cnt,cnt+nklist] = -Layer.t_orth[ilay]
                    elseif ilay==Layer.nlayer
                        H[cnt,cnt-(cnt_perlay[ilay]-cnt_perlay[ilay-1])] = -Layer.t_orth[ilay-1]
                    else
                        H[cnt,cnt+nklist] = -Layer.t_orth[ilay]
                        H[cnt,cnt-(cnt_perlay[ilay]-cnt_perlay[ilay-1])] = -Layer.t_orth[ilay-1]
                    end
                end
                cnt += 1
            end
        end
    end
    return H
end



"""
    getRhoatmu(Layer::LayerParameters,norb::Int64,orbs_layer::Vector{Int64},k::Vector{Float64},
                klist::Vector{Vector{Float64}},eta::Float64)

Compute the spectral weight at the fermi level at k point from the layers' parameters in Layer. The
size of Hamiltonian is given by norb, and there is (2*)ncdw k-points per layer.
"""
function getRhoatmu(Layer::LayerParameters,norb::Int64,orbs_layer::Vector{Int64},k::Vector{Float64},
                    klist::Vector{Vector{Float64}}, eta::Float64)

    rho_out = zeros(Float64,norb)
    #Shift klist
    klistH = [[mod((klist[ik]+k)[1]+pi,2pi)-pi,mod((klist[ik]+k)[2]+pi,2pi)-pi] for ik=1:length(klist)]
    #Get Hamiltonian
    H = getH(Layer,norb,klistH)
    #evaluate Green's function at mu
    rho = -imag(inv(1im*eta*I(norb)-H))/pi
    #Extract diagonal and eventually scale by gt
    for iorb=1:norb
        rho_out[iorb]=gt(Layer.x[orbs_layer[iorb]])*rho[iorb,iorb]
    end
    return rho_out
end

"""
    getBareFS(Layer::LayerParam,nk::Int64,eta::Float64)

Calculates the Fermi Surface without magnetic field from the Layers' parameters
Layer, the number of k-points nk.
"""
function getBareFS(Layer::LayerParameters,nk::Int64,eta::Float64)
    #Initialize FS struct
    norb,orbs,orbs_layer,klist = prepFSHam(Layer)
    FS = FermiSurface(nk=nk,norb=norb,orbs=orbs)
    #Compute FS at each k-point
    for ikx = 1:nk
        for iky = 1:nk
            FS.fs[ikx,iky,:] = getRhoatmu(Layer,norb,orbs_layer,FS.kgrid[ikx,iky,:],klist,eta)
        end
    end
    return FS
end
