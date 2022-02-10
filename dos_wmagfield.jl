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

    dos         :: FieldMesh

    DensityOfStates(;dos_size) = new(zeros(Float64,dos_size))

end




################# Functions #################################################

"""
    getBdepDOS(pSim::SimulationParameters,Layer::LayerParameters)

main routine that manage the DOS calculation over the magnetic field grid

Magnetic field is rotated by pSim.theta angle, and we put ourselves in the Landau gauge chosen
so that we can perform FT along the y direction.

"""
function getBdepDOS(pSim::SimulationParameters,Layer::LayerParameters)
    c = prepareHparam(Layer)
end

"""
    prepareHparam(args)

Generates the essential parameters for the blocks of the Hamiltonian : size, orbitals, etc..
"""
function prepareHparam(Layer::LayerParameters)
    #First search for CDW
    ncdw, ncdw_list, qlist = get_ncdw_wMagField(Layer)
    #Get the k-vector list
    klist = get_klist_wMagField(ncdw,ncdw_list,qlist)

    #Compute norb by checking if YRZ & construct orbs array
    norb=0
    orbs_layer = Int64[]
    orb_cnt=1
    for ilayer = 1:nlayer
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
    for ilayer = 1:nlayer
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
        return 0,[],[]
    end
    #For all different Q, find the necessary nbr of kpoints -> only for Qy since no FT along x
    ncdw_list = []
    for q in qlist
        if abs(q[2])>1e-8 && abs(q[1])<1e-8           #if y unidirectional
            append!(ncdw_list,Int64(2*pi/abs(q[2])))
        elseif abs(q[1])>1e-8 && abs(q[2])>1e-8           #if (x,y) bi-directional
            throw(ArgumentError("Q=(qx!=0,qy!=0) not yet implemented for QO calculations."))
        end
    end
    #Total nbr is the product
    return Int64(prod(ncdw_list)),ncdw_list,qlist
end

"""
    get_klist(ncdw::Int64,ncdw_list::Vector{Any},qlist::Vector{Any})

Determines the list of k-points necessary for FS calculation from the q vectors
"""
function get_klist_wMagField(ncdw::Int64,ncdw_list::Vector{Any},qlist::Vector{Any})
    klist=[[0.0,0.0]]
    #If no cdw vector
    if ncdw == 0
        return klist
    end
    #Else, find all possible k-points
    ktmp = [[0.0,0.0]]
    for iq=1:length(qlist)
        for k in ktmp
            for inq=1:ncdw_list[iq]
                if !([mod((k+qlist[iq])[1],2pi),mod((k+qlist[iq])[2],2pi)] in klist)
                    append!(klist,[k+qlist[iq]])
                end
            end
        end
        if length(ktmp)==length(klist)
            if iq==length(qlist)
                break
            else
                @goto loopagain
            end
        else
            ktmp=klist
        end
        @label loopagain
    end
    return klist
end
