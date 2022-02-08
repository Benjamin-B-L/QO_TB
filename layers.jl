###Type defs
const LayerParam = Array{Float64,1}
const Qvector = Dict{Int64, Vector}


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
- nlayer          : Number of layers considered
- t_orth          : interlayer coupling (if nlayer>1)
- ac              : interlayer distance (in units of a=1)

end
"""

struct LayerParameters
    nlayer  :: Int64
    t0      :: LayerParam
    tp0     :: LayerParam
    t       :: LayerParam
    tp      :: LayerParam
    J       :: LayerParam
    chi     :: LayerParam
    x       :: LayerParam
    D       :: LayerParam
    mu      :: LayerParam
    mu_aux  :: LayerParam
    Q       :: Qvector
    t_orth  :: LayerParam
    ac      :: LayerParam

    function LayerParameters(;nlayer,t0,tp0,J=zeros(Float64,nlayer),chi=zeros(Float64,nlayer),
                                x=ones(Float64,nlayer),D=zeros(Float64,nlayer),mu=zeros(Float64,nlayer),
                                mu_aux=zeros(Float64,nlayer),Q=Dict(i=>[] for i=1:nlayer),
                                t_orth=Array{Float64}[],ac=Array{Float64}[])

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
        if length(mu) != nlayer
            push!(warn,"mu")
        end
        if length(mu_aux) != nlayer
            push!(warn,"mu_aux")
        end
        if length(Q) != nlayer
            push!(warn,"Q")
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

        ## generate struct
        return new(nlayer,t0,tp0,t0,tp0,J,chi,x,D,mu,mu_aux,Q,t_orth,ac)
    end
end
