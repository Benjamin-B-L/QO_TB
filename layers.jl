###Type defs
const LayerParam = Array{Float64,1}
const Qvector = Dict{Int64, Vector}
const ModeParam = Array{String,1}


########## Layers struct ###########

"""
    LayerParameters
    Contains the parameters for each layer of the model

Fields:    all of them but nlayer, t_orth and ac are arrays of length nlayer
----------

- t0              : Bare NN hopping
- tp0             : Bare NNN hopping
- t               : Renormalized NN hopping
- tp              : Renormalized NNN hopping
- J               : Exchange parameter for hopping renormalization
- chi             : Susceptibility for hopping renormalization
- x               : doping for hopping renormalization
- D               : Coupling between physical and auxiliary electrons
- mu              : chemical potential for the physical electrons
- mu_aux          : chemical potential for the auxiliary electrons
- Q               : CDW wave vector (set to 0 if no CDW)
- Pcdw            : CDW coupling constant
- nlayer          : Number of layers considered
- t_orth          : interlayer coupling (if nlayer>1)
- ac              : interlayer distance (in units of a=1)

end
"""

struct LayerParameters
    nlayer   :: Int64
    t0       :: LayerParam
    tp0      :: LayerParam
    t        :: LayerParam
    tp       :: LayerParam
    J        :: LayerParam
    chi      :: LayerParam
    x        :: LayerParam
    D        :: LayerParam
    yrz_mode :: ModeParam
    mu       :: LayerParam
    mu_aux   :: LayerParam
    Q        :: Qvector
    Pcdw     :: LayerParam
    t_orth   :: LayerParam
    ac       :: LayerParam

    function LayerParameters(;nlayer,t0,tp0,J=zeros(Float64,nlayer),chi=zeros(Float64,nlayer),
                                x=ones(Float64,nlayer),D=zeros(Float64,nlayer),yrz_mode=["Norm" for i=1:nlayer]
                                ,mu=zeros(Float64,nlayer),mu_aux=zeros(Float64,nlayer),Q=Dict(i=>[] for i=1:nlayer),
                                Pcdw=zeros(Float64,nlayer),t_orth=zeros(Float64,nlayer-1),ac=zeros(Float64,nlayer-1))

        ## Check that all variables are consistent with nlayer
        warn=[]
        if length(t0) != nlayer
            push!(warn,"t0")
        end
        if length(tp0) != nlayer
            push!(warn,"tp0")
        end
        if length(J) != nlayer
            push!(warn,"J")
        end
        if length(chi) != nlayer
            push!(warn,"chi")
        end
        if length(x) != nlayer
            push!(warn,"x")
        end
        if length(D) != nlayer
            push!(warn,"D")
        end
        if length(yrz_mode) != nlayer
            push!(warn,"yrz_mode")
        end
        if length(mu) != nlayer
            push!(warn,"mu")
        end
        if length(mu_aux) != nlayer
            push!(warn,"mu_aux")
        end
        if length(Q) != nlayer
            push!(warn,"Q")
        end
        if length(Pcdw) != nlayer
            push!(warn,"Pcdw")
        end
        if length(t_orth) != nlayer-1
            push!(warn,"t_orth")
        end
        if length(ac) != nlayer-1
            push!(warn,"ac")
        end
        if length(warn)>0
            throw(ArgumentError("The following arrays should be of the same size as nlayer : $warn"))
        end

        #Pre-define renormalized hoppings
        t = [t0[ilay] for ilay=1:nlayer]
        tp = [tp0[ilay] for ilay=1:nlayer]

        for ilay=1:nlayer
            #Check that CDW coupling is zero if no Q vector
            if abs(Pcdw[ilay])>1e-8 && length(Q[ilay])==0
                throw(ArgumentError("If no Q vector is provided for a layer, Pcdw should be set to zero !"))
            end
            #renormalize hopping
            t[ilay] = t0[ilay]*gt(x[ilay]) + 3*gj(x[ilay])*J[ilay]*chi[ilay]/8
            tp[ilay] = tp0[ilay]*gt(x[ilay])
        end

        ## generate struct
        return new(nlayer,t0,tp0,t,tp,J,chi,x,D,yrz_mode,mu,mu_aux,Q,Pcdw,t_orth,ac)
    end
end




################# Functions #################################################

"""
    gt(x)

calculates function gt entering renormalization factor of hopping terms
"""
function gt(x)
    return (2*x)/(1+x)
end

"""
    gj(x)

calculates function gj entering renormalization factor of hopping terms
"""
function gj(x)
    return 4/((1+x)^2)
end
