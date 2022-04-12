using LinearAlgebra

include("params.jl")
include("layers.jl")


########## DOS struct ###########

"""
    DensityOfStates
    Contains the calculated density of states on the magnetic field grid.
    TO DO : will also be used to store QO spectrum

Fields:
----------

- dos             : density of states
"""

struct DensityOfStates

    dos         :: Matrix{Float64}

    DensityOfStates(;dos_size,nky) = new(zeros(Float64,dos_size,nky))

end




################# Functions #################################################

"""
    getBdepDOS(pSim::SimulationParameters,Layer::LayerParameters)

main routine that manage the DOS calculation over the magnetic field grid

Magnetic field is rotated by pSim.theta angle, and we put ourselves in the Landau gauge chosen
so that we can perform FT along the y direction.

"""
function getBdepDOS(pSim::SimulationParameters,Layer::LayerParameters)
    #Prepare parameters for the Hamiltonian blocks
    norb,orbs_layer,klist = prepareHparam(Layer)

    dos = DensityOfStates(dos_size=pSim.nb,nky=pSim.ny_avg+1)
    for iy = 1:pSim.ny_avg
        #Pick a random ky between [0:2pi[ and create ky_list
        ky = 2*pi*rand(Float64)
        ky_list = [mod(klist[ik][2]+ky,2*pi) for ik=1:length(klist)]
        #Loop over magnetic field values
        for ib = 1:pSim.nb
            diag,diagL,diagR = getDiags(pSim.nx,pSim.bgrid[ib],pSim.theta,ky_list,Layer,norb,0.5,pSim.eta)
            dos.dos[ib,iy] = get_dos(diag,diagL,diagR,Layer,pSim.nx,length(ky_list))
        end
    end
    for ib=1:pSim.nb
        dos.dos[ib,pSim.ny_avg+1] = sum(dos.dos[ib,:])/pSim.ny_avg
    end
    return dos
end

"""
    prepareHparam(args)

Generates the essential parameters for the blocks of the Hamiltonian : size, orbitals, etc..
"""
function prepareHparam(Layer::LayerParameters)
    #First search for CDW
    ncdw_list, qlist = get_ncdw_wMagField(Layer)
    #Get the k-vector list
    klist = get_klist_wMagField(ncdw_list,qlist)
    ncdw = length(klist)

    #Compute norb by checking if YRZ & construct orbs array
    norb=0
    orbs_layer = Int64[]
    orb_cnt=1
    for ilayer = 1:Layer.nlayer
        if Layer.D[ilayer] < 1e-8
            norb += max(ncdw,1)              #No auxiliary fermions in layer ilayer
            for icdw = 1:max(ncdw,1)         #get all orbs
                append!(orbs_layer,ilayer)
            end
        else
            norb += 2*max(ncdw,1)            #Auxiliary fermions in layer ilayer
            for icdw = 1:max(ncdw,1)
                append!(orbs_layer,ilayer)
            end
            for icdw = 1:max(ncdw,1)
                append!(orbs_layer,ilayer)
            end
        end
    end

    return norb,orbs_layer,klist
end

"""
    get_ncdw(Layer::LayerParameters)

Determines the number of different k-points necessary from all Q vectors in all layers
Note that contrary to bare FS calculation, we keep only the (0,Qy) vectors because FT is only
done in y direction.
"""
function get_ncdw_wMagField(Layer::LayerParameters)
    #Extract all different Q vectors
    qlist = []
    for ilayer = 1:Layer.nlayer
        if length(Layer.Q[ilayer]) > 0
            for iq = 1:length(Layer.Q[ilayer])
                if !(Layer.Q[ilayer][iq] in qlist)
                    if abs(Layer.Q[ilayer][iq][1])<1e-8 && abs(Layer.Q[ilayer][iq][2])>1e-8
                        append!(qlist,[Layer.Q[ilayer][iq]])
                    end
                end
            end
        end
    end

    #If no cdw vector
    if length(qlist)==0
        return [],[],[]
    end
    #For all different Q, find the necessary nbr of kpoints -> only for Qy since no FT along x
    ncdw_list = []
    for q in qlist
        if abs(q[2])>1e-8 && abs(q[1])<1e-8           #if y unidirectional
            if abs(trunc(Int64,2*pi/abs(q[2]))-2*pi/abs(q[2]))<1e-4
                append!(ncdw_list,trunc(Int64,2*pi/abs(q[2])))
            else
                throw(ArgumentError("Q vector provided is not an integer..."))
            end
        elseif abs(q[1])>1e-8 && abs(q[2])>1e-8           #if (x,y) bi-directional
            throw(ArgumentError("Q=(qx!=0,qy!=0) not yet implemented for QO calculations."))
        end
    end
    #Total nbr is the product
    return ncdw_list,qlist
end

"""
    get_klist(ncdw::Int64,ncdw_list::Vector{Any},qlist::Vector{Any})

Determines the list of k-points necessary for FS calculation from the q vectors
"""
function get_klist_wMagField(ncdw_list::Vector{Any},qlist::Vector{Any})
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
    getDiags(nx::Int64,B::Float64,theta::Float64,ky_list::Vector{Float64},Layer::LayerParameters,norb::Int64,sigma::Float64,eta::Float64)

Compute iteratively the left and right diagonal blocks necessary for the Green's function calculation
"""
function getDiags(nx::Int64,B::Float64,theta::Float64,ky_list::Vector{Float64},Layer::LayerParameters,norb::Int64,sigma::Float64,eta::Float64)
    diag = zeros(Complex{Float64},nx,norb,norb)
    diagL = zeros(Complex{Float64},nx,norb,norb)
    diagR = zeros(Complex{Float64},nx,norb,norb)

    #Diag
    for ix = 1:nx
        diag[ix,:,:] = diagblock(ix,B,theta,ky_list,Layer,norb,sigma,eta)
    end
    #Left Diag
    diagL[1,:,:] = diagblock(1,B,theta,ky_list,Layer,norb,sigma,eta)
    #Right Diag
    diagR[nx,:,:] = diagblock(nx,B,theta,ky_list,Layer,norb,sigma,eta)
    for ix=2:nx
        #Left Diag
        diagL[ix,:,:] = diag[ix,:,:]-offdiagblock(ix,ix-1,B,theta,ky_list,Layer,norb,sigma)*inv(diagL[ix-1,:,:])*offdiagblock(ix-1,ix,B,theta,ky_list,Layer,norb,sigma)
        #Right Diag
        diagR[nx-ix+1,:,:] = diag[nx-ix+1,:,:]-offdiagblock(nx-ix+1,nx-ix+1+1,B,theta,ky_list,Layer,norb,sigma)*inv(diagR[nx-ix+1+1,:,:])*offdiagblock(nx-ix+1+1,nx-ix+1,B,theta,ky_list,Layer,norb,sigma)
    end

    return diag,diagL,diagR
end

"""
    diagblock(ix::Int64,B::Float64,theta::Float64,ky_list::Vector{Float64},Layer::LayerParameters,
               norb::Int64,sigma::Float64)

Calculate the diagonal block n° ix
"""
function diagblock(ix::Int64,B::Float64,theta::Float64,ky_list::Vector{Float64},Layer::LayerParameters,
                    norb::Int64,sigma::Float64,eta::Float64)

    nklist = length(ky_list)
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
                        #Dispersion
                        H[cnt,cnt] += eps("diag",ky_list[ik],Layer.t[ilay],0.0,Layer.mu[ilay],ix,0.0,B,theta)
                        #YRZ coupling
                        H[cnt,cnt+nklist] += yrz_delta("diag",ky_list[ik],Layer.D[ilay],ix,B,theta)
                        #CDW
                        for q in Layer.Q[ilay]
                            if abs(q[1])<1e-8
                                #+Qy
                                kq = mod((ky_list[ik]+q[2]),2pi)
                                kq_index = cnt_perlay[ilay] + findfirst(x->abs(x-kq)^2<1e-8,ky_list)
                                H[cnt,kq_index] += cdw_coupling("diag","y",ky_list[ik],Layer.Pcdw[ilay],ix,0.0,B,theta)
                                H[kq_index,cnt] += cdw_coupling("diag","y",ky_list[ik],Layer.Pcdw[ilay],ix,0.0,B,theta)
                                #-Qy
                                kq = mod((ky_list[ik]-q[2]),2pi)
                                kq_index = cnt_perlay[ilay] + findfirst(x->abs(x-kq)^2<1e-8,ky_list)
                                H[cnt,kq_index] += cdw_coupling("diag","y",ky_list[ik],Layer.Pcdw[ilay],ix,0.0,B,theta)
                                H[kq_index,cnt] += cdw_coupling("diag","y",ky_list[ik],Layer.Pcdw[ilay],ix,0.0,B,theta)
                            else
                                #+Qx
                                H[cnt,cnt] += cdw_coupling("diag","x",ky_list[ik],Layer.Pcdw[ilay],ix,0.0,B,theta,q[1])
                                #-Qx
                                H[cnt,cnt] += cdw_coupling("diag","x",ky_list[ik],Layer.Pcdw[ilay],ix,0.0,B,theta,-q[1])
                            end
                        end
                        #interlayer coupling
                        if Layer.nlayer > 1
                            if ilay==1
                                H[cnt,cnt+2*nklist] = -Layer.t_orth[ilay]*exp(im*2*pi*ix*B*cos(theta)*tan(Layer.ac[ilay]))
                            elseif ilay==Layer.nlayer
                                H[cnt,cnt-(cnt_perlay[ilay]-cnt_perlay[ilay-1])] = -Layer.t_orth[ilay-1]*exp(-im*2*pi*ix*B*cos(theta)*tan(Layer.ac[ilay-1]))
                            else
                                H[cnt,cnt+2*nklist] = -Layer.t_orth[ilay]*exp(im*2*pi*ix*B*cos(theta)*tan(Layer.ac[ilay]))
                                H[cnt,cnt-(cnt_perlay[ilay]-cnt_perlay[ilay-1])] = -Layer.t_orth[ilay-1]*exp(-im*2*pi*ix*B*cos(theta)*tan(Layer.ac[ilay-1]))
                            end
                        end
                        #Zeeman splitting
                        H[cnt,cnt] += 4*pi*Layer.g_tilde[ilay]*B*sigma
                    #Auxiliary electrons
                    else
                        #Dispersion
                        H[cnt,cnt] += eps_aux("diag",ky_list[ik],Layer.t[ilay],0.0,Layer.mu_aux[ilay],Layer.yrz_mode[ilay],ix,0.0,B,theta)
                        #YRZ
                        H[cnt,cnt-nklist] += conj(yrz_delta("diag",ky_list[ik],Layer.D[ilay],ix,B,theta))
                        #Zeeman splitting
                        H[cnt,cnt] += 4*pi*Layer.g_tilde[ilay]*B*sigma
                    end
                    cnt += 1
                end
            end
        #if no YRZ
        else
            #Only physical electrons
            for ik = 1:nklist
                #Diagonal
                H[cnt,cnt] += eps("diag",ky_list[ik],Layer.t[ilay],0.0,Layer.mu[ilay],ix,0.0,B,theta)
                #CDW
                for q in Layer.Q[ilay]
                    if abs(q[1])<1e-8
                        #+Qy
                        kq = mod((ky_list[ik]+q[2]),2pi)
                        kq_index = cnt_perlay[ilay] + findfirst(x->abs(x-kq)^2<1e-8,ky_list)
                        H[cnt,kq_index] += cdw_coupling("diag","y",ky_list[ik],Layer.Pcdw[ilay],ix,0.0,B,theta)
                        H[kq_index,cnt] += cdw_coupling("diag","y",ky_list[ik],Layer.Pcdw[ilay],ix,0.0,B,theta)
                        #-Qy
                        kq = mod((ky_list[ik]-q[2]),2pi)
                        kq_index = cnt_perlay[ilay] + findfirst(x->abs(x-kq)^2<1e-8,ky_list)
                        H[cnt,kq_index] += cdw_coupling("diag","y",ky_list[ik],Layer.Pcdw[ilay],ix,0.0,B,theta)
                        H[kq_index,cnt] += cdw_coupling("diag","y",ky_list[ik],Layer.Pcdw[ilay],ix,0.0,B,theta)
                    else
                        #+Qx
                        H[cnt,cnt] += cdw_coupling("diag","x",ky_list[ik],Layer.Pcdw[ilay],ix,0.0,B,theta,q[1])
                        #-Qx
                        H[cnt,cnt] += cdw_coupling("diag","x",ky_list[ik],Layer.Pcdw[ilay],ix,0.0,B,theta,-q[1])
                    end
                end
                #interlayer coupling
                if Layer.nlayer > 1
                    if ilay==1
                        H[cnt,cnt+nklist] = -Layer.t_orth[ilay]*exp(im*2*pi*ix*B*cos(theta)*tan(Layer.ac[ilay]))
                    elseif ilay==Layer.nlayer
                        H[cnt,cnt-(cnt_perlay[ilay]-cnt_perlay[ilay-1])] = -Layer.t_orth[ilay-1]*exp(-im*2*pi*ix*B*cos(theta)*tan(Layer.ac[ilay-1]))
                    else
                        H[cnt,cnt+nklist] = -Layer.t_orth[ilay]*exp(im*2*pi*ix*B*cos(theta)*tan(Layer.ac[ilay]))
                        H[cnt,cnt-(cnt_perlay[ilay]-cnt_perlay[ilay-1])] = -Layer.t_orth[ilay-1]*exp(-im*2*pi*ix*B*cos(theta)*tan(Layer.ac[ilay-1]))
                    end
                end
                #Zeeman splitting
                H[cnt,cnt] += 4*pi*Layer.g_tilde[ilay]*B*sigma
                cnt += 1
            end
        end
    end
    return im*eta*I(norb)-H
end


"""
    offdiagblock(ix1::Int64,ix2::Int64,B::Float64,theta::Float64,ky_list::Vector{Float64},Layer::LayerParameters,
               norb::Int64,sigma::Float64)

Calculate the offdiagonal block n° ix
"""
function offdiagblock(ix1::Int64,ix2::Int64,B::Float64,theta::Float64,ky_list::Vector{Float64},Layer::LayerParameters,
                    norb::Int64,sigma::Float64)

    nklist = length(ky_list)
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
                        #Dispersion
                        H[cnt,cnt] += eps("offdiag",ky_list[ik],Layer.t[ilay],Layer.tp[ilay],Layer.mu[ilay],ix1,ix2,B,theta)
                        #YRZ coupling
                        H[cnt,cnt+nklist] += yrz_delta("offdiag",ky_list[ik],Layer.D[ilay],0,B,theta)
                        #CDW
                        for q in Layer.Q[ilay]
                            if abs(q[1])<1e-8
                                #+Qy
                                kq = mod((ky_list[ik]+q[2]),2pi)
                                kq_index = cnt_perlay[ilay] + findfirst(x->abs(x-kq)^2<1e-8,ky_list)
                                H[cnt,kq_index] += cdw_coupling("offdiag","y",ky_list[ik],Layer.Pcdw[ilay],ix1,0.0,B,theta)
                                H[kq_index,cnt] += cdw_coupling("offdiag","y",ky_list[ik],Layer.Pcdw[ilay],ix1,0.0,B,theta)
                                #-Qy
                                kq = mod((ky_list[ik]-q[2]),2pi)
                                kq_index = cnt_perlay[ilay] + findfirst(x->abs(x-kq)^2<1e-8,ky_list)
                                H[cnt,kq_index] += cdw_coupling("offdiag","y",ky_list[ik],Layer.Pcdw[ilay],ix1,0.0,B,theta)
                                H[kq_index,cnt] += cdw_coupling("offdiag","y",ky_list[ik],Layer.Pcdw[ilay],ix1,0.0,B,theta)
                            else
                                #+Qx
                                H[cnt,cnt] += cdw_coupling("offdiag","x",ky_list[ik],Layer.Pcdw[ilay],ix1,ix2,B,theta,q[1])
                                #-Qx
                                H[cnt,cnt] += cdw_coupling("offdiag","x",ky_list[ik],Layer.Pcdw[ilay],ix1,ix2,B,theta,-q[1])
                            end
                        end
                    #Auxiliary electrons
                    else
                        #Dispersion
                        H[cnt,cnt] += eps_aux("offdiag",ky_list[ik],Layer.t[ilay],Layer.tp[ilay],Layer.mu_aux[ilay],Layer.yrz_mode[ilay],ix1,ix2,B,theta)
                        #YRZ
                        H[cnt,cnt-nklist] += conj(yrz_delta("offdiag",ky_list[ik],Layer.D[ilay],0,B,theta))
                    end
                    cnt += 1
                end
            end
        #if no YRZ
        else
            #Only physical electrons
            for ik = 1:nklist
                #Diagonal
                H[cnt,cnt] += eps("offdiag",ky_list[ik],Layer.t[ilay],Layer.tp[ilay],Layer.mu[ilay],ix1,ix2,B,theta)
                #CDW
                for q in Layer.Q[ilay]
                    if abs(q[1])<1e-8
                        #+Qy
                        kq = mod((ky_list[ik]+q[2]),2pi)
                        kq_index = cnt_perlay[ilay] + findfirst(x->abs(x-kq)<1e-8,ky_list)
                        H[cnt,kq_index] += cdw_coupling("offdiag","y",ky_list[ik],Layer.Pcdw[ilay],ix1,0.0,B,theta)
                        H[kq_index,cnt] += cdw_coupling("offdiag","y",ky_list[ik],Layer.Pcdw[ilay],ix1,0.0,B,theta)
                        #-Qy
                        kq = mod((ky_list[ik]-q[2]),2pi)
                        kq_index = cnt_perlay[ilay] + findfirst(x->abs(x-kq)^2<1e-8,ky_list)
                        H[cnt,kq_index] += cdw_coupling("offdiag","y",ky_list[ik],Layer.Pcdw[ilay],ix1,0.0,B,theta)
                        H[kq_index,cnt] += cdw_coupling("offdiag","y",ky_list[ik],Layer.Pcdw[ilay],ix1,0.0,B,theta)
                    else
                        #+Qx
                        H[cnt,cnt] += cdw_coupling("offdiag","x",ky_list[ik],Layer.Pcdw[ilay],ix1,ix2,B,theta,q[1])
                        #-Qx
                        H[cnt,cnt] += cdw_coupling("offdiag","x",ky_list[ik],Layer.Pcdw[ilay],ix1,ix2,B,theta,-q[1])
                    end
                end
                cnt += 1
            end
        end
    end
    return -H
end


"""
    eps(switch,ky,t,tp,ix,B,theta)

Compute the dispersion term for the diagonal or off-diag block for physical electrons.
"""
function eps(switch,ky,t,tp,mu,ix1,ix2,B,theta)
    if switch=="diag"
        return -2*t*cos(2*pi*B*cos(theta)*ix1+ky)-mu
    elseif switch=="offdiag"
        return -t-2*tp*cos(2*pi*B*cos(theta)*(ix1+0.5*sign(ix2-ix1))+ky)
    else
        throw(ArgumentError("Input switch not recongnized: must be either set to diag or offdiag"))
    end
end

"""
    eps_aux(switch,ky,t,tp,ix,B,theta)

Compute the dispersion term for the diagonal or off-diag block for auxiliary electrons.
"""
function eps_aux(switch,ky,t,tp,mu,yrz_mode,ix1,ix2,B,theta)
    if switch=="diag"
        return 2*t*cos(2*pi*B*cos(theta)*ix1+ky)-mu
    elseif switch=="offdiag"
        if yrz_mode=="Norm"
            return t
        elseif yrz_mode=="AFM"
            return t-2*tp*cos(2*pi*B*cos(theta)*(ix1+0.5*sign(ix2-ix1))+ky)
        else
            throw(ArgumentError("YRZ mode $yrz_mode not implemented"))
        end
    else
        throw(ArgumentError("Input switch not recongnized: must be either set to diag or offdiag"))
    end
end


"""
    yrz_delta(switch,ky,t,ix,B,theta)

Compute the YRZ coupling term.
"""
function yrz_delta(switch,ky,D,ix,B,theta)
    if switch=="diag"
        return -D*cos(2*pi*B*cos(theta)*ix+ky)
    elseif switch=="offdiag"
        return 0.5*D
    else
        throw(ArgumentError("Input switch not recongnized: must be either set to diag or offdiag"))
    end
end

"""
    cdw_coupling(switch,direction,ky,Pcdw,ix,B,theta)

Compute the CDW coupling.
"""
function cdw_coupling(switch,direction,ky,Pcdw,ix1,ix2,B,theta,qx=0)
    if switch=="diag"
        if direction=="y"
            return -Pcdw*cos(2*pi*B*cos(theta)*ix1+ky)
        elseif direction=="x"
            return -2*Pcdw*cos(2*pi*B*cos(theta)*ix1+ky)*cos(qx*ix1)
        end
    elseif switch=="offdiag"
        if direction=="y"
            return 0.5*Pcdw
        elseif direction=="x"
            return 0.5*Pcdw*(exp(im*qx*(ix1+sign(ix2-ix1)))+exp(-im*qx*ix1))
        end
    else
        throw(ArgumentError("Input switch not recongnized: must be either set to diag or offdiag"))
    end
end

"""
    get_dos(diag::Array{ComplexF64, 3},diagL::Array{ComplexF64, 3},diagR::Array{ComplexF64, 3},
    Layer::LayerParameters,nx::Int64,nklist::Int64)

documentation
"""
function get_dos(diag::Array{ComplexF64, 3},diagL::Array{ComplexF64, 3},diagR::Array{ComplexF64, 3},
                  Layer::LayerParameters,nx::Int64,nklist::Int64)
    dos = 0

    for ix=1:nx
        cnt = 1
        G = inv(-diag[ix,:,:]+diagL[ix,:,:]+diagR[ix,:,:])
        for ilay=1:Layer.nlayer
            #Sum over physical electrons
            for ik=1:nklist
                dos += gt(Layer.x[ilay])*G[cnt,cnt]
                cnt +=1
            end
            #Pass over auxiliary electrons if any
            if abs(Layer.D[ilay])>1e-8
                cnt+=nklist
            end
        end
    end
    return -imag(dos)/(pi*nx)

end
