include("params.jl")
include("layers.jl")

######### INITIALIZATION ##########

"""
    initialize(file)

Generates the SimulationParameters and LayerParameters structs from the input file.
"""
function initialize(file::String)
    ###Check if file exist
    if isfile(file)==false
        throw(ArgumentError("Cannot find input file $file !"))
    end
    input = readlines(file)
    #Remove comments
    for iline=1:length(input)
        input[iline] = split(input[iline],"#")[1]
    end
    ###SimulationParameters
    #nx
    indx = findfirst(x->contains(x,"nx"),input)
    if indx==nothing
        throw(ArgumentError("missing parameter nx !"))
    end
    nx = parse(Int64,split(input[indx],"=")[2])
    #ny_avg
    indx = findfirst(x->contains(x,"ny_avg"),input)
    if indx==nothing
        throw(ArgumentError("missing parameter ny_avg !"))
    end
    ny_avg = parse(Int64,split(input[indx],"=")[2])
    #bmin
    indx = findfirst(x->contains(x,"bmin"),input)
    if indx==nothing
        throw(ArgumentError("missing parameter bmin !"))
    end
    tmp = split(input[indx],"=")[2]
    if contains(tmp,"/")
        bmin = parse(Float64,split(tmp,"/")[1])/parse(Float64,split(tmp,"/")[2])
    else
        bmin = parse(Float64,split(input[indx],"=")[2])
    end
    #bmax
    indx = findfirst(x->contains(x,"bmax"),input)
    if indx==nothing
        throw(ArgumentError("missing parameter bmax !"))
    end
    tmp = split(input[indx],"=")[2]
    if contains(tmp,"/")
        bmax = parse(Float64,split(tmp,"/")[1])/parse(Float64,split(tmp,"/")[2])
    else
        bmax = parse(Float64,tmp)
    end
    #nb
    indx = findfirst(x->contains(x,"nb"),input)
    if indx==nothing
        throw(ArgumentError("missing parameter nb !"))
    end
    nb = parse(Int64,split(input[indx],"=")[2])
    #theta_in
    indx = findfirst(x->contains(x,"theta_in"),input)
    if indx==nothing
        throw(ArgumentError("missing parameter theta_in !"))
    end
    theta_in = parse(Float64,split(input[indx],"=")[2])
    #eta
    indx = findfirst(x->contains(x,"broadening"),input)
    if indx==nothing
        throw(ArgumentError("missing parameter broadening!"))
    end
    eta = parse(Float64,split(input[indx],"=")[2])
    #nk
    indx = findfirst(x->contains(x,"nk"),input)
    if indx==nothing
        throw(ArgumentError("missing parameter nk !"))
    end
    nk = parse(Int64,split(input[indx],"=")[2])
    #outfolder
    indx = findfirst(x->contains(x,"outfolder"),input)
    if indx==nothing
        throw(ArgumentError("missing parameter outfolder !"))
    end
    outfolder = String(split(input[indx],"=")[2])

    pSim = SimulationParameters(nx=nx,ny_avg=ny_avg,bmin=bmin,bmax=bmax,nb=nb,theta_in=theta_in,eta=eta,nk=nk,outfolder=outfolder)

    ###LayerParameters
    missing_input=[]
    #nlayer
    indx = findfirst(x->contains(x,"nlayer"),input)
    if indx==nothing
        throw(ArgumentError("missing parameter nlayer !"))
    end
    nlayer = parse(Int64,split(input[indx],"=")[2])
    #t0
    indx = findfirst(x->contains(x,"t0"),input)
    if indx==nothing
        throw(ArgumentError("missing parameter t0 !"))
    end
    t0 = [parse(Float64,x) for x in split(split(input[indx],"=")[2])]
    #tp0
    indx = findfirst(x->contains(x,"tp0"),input)
    if indx==nothing
        throw(ArgumentError("missing parameter tp0 !"))
    end
    tp0 = [parse(Float64,x) for x in split(split(input[indx],"=")[2])]
    #J
    indx = findfirst(x->contains(x,"J"),input)
    if indx==nothing
        push!(missing_input,"J")
        J = zeros(Float64,nlayer)
    else
        J = [parse(Float64,x) for x in split(split(input[indx],"=")[2])]
    end
    #chi
    indx = findfirst(x->contains(x,"chi"),input)
    if indx==nothing
        push!(missing_input,"chi")
        chi = zeros(Float64,nlayer)
    else
        chi = [parse(Float64,x) for x in split(split(input[indx],"=")[2])]
    end
    #x
    indx = findfirst(x->contains(x,"doping"),input)
    if indx==nothing
        push!(missing_input,"doping")
        x = ones(Float64,nlayer)
    else
        x = [parse(Float64,x) for x in split(split(input[indx],"=")[2])]
    end
    #D
    indx = findfirst(x->contains(x,"D"),input)
    if indx==nothing
        push!(missing_input,"D")
        D = zeros(Float64,nlayer)
    else
        D = [parse(Float64,x) for x in split(split(input[indx],"=")[2])]
    end
    #yrz_mode
    indx = findfirst(x->contains(x,"yrz_mode"),input)
    if indx==nothing
        push!(missing_input,"yrz_mode")
        yrz_mode = ["Norm" for ilay=1:nlayer]
    else
        yrz_mode = [String(x) for x in split(split(input[indx],"=")[2])]
    end
    #mu_phys
    indx = findfirst(x->contains(x,"mu_phys"),input)
    if indx==nothing
        push!(missing_input,"mu_phys")
        mu = zeros(Float64,nlayer)
    else
        mu = [parse(Float64,x) for x in split(split(input[indx],"=")[2])]
    end
    #mu_aux
    indx = findfirst(x->contains(x,"mu_aux"),input)
    if indx==nothing
        push!(missing_input,"mu_aux")
        mu_aux = zeros(Float64,nlayer)
    else
        mu_aux = [parse(Float64,x) for x in split(split(input[indx],"=")[2])]
    end
    #Q
    indx = findfirst(x->contains(x,"Q"),input)
    if indx==nothing
        push!(missing_input,"Q")
        Q=Dict(i=>[] for i=1:nlayer)
    else
        Q = Dict{Int64,Vector}()
        cnt = 1
        for x in split(split(input[indx],"=")[2])
            vec = []
            for q in SubString.(x,findall(r"\(.*?\)",x))
                push!(vec,[parse(Float64,SubString.(q,findall(r"\d\.\d{1,8}",q))[1]),parse(Float64,SubString.(q,findall(r"\d\.\d{1,8}",q))[2])])
            end
            merge!(Q,Dict{Int64,Vector}(cnt=>vec))
            cnt+=1
        end
    end
    #Pcdw
    indx = findfirst(x->contains(x,"Pcdw"),input)
    if indx==nothing
        push!(missing_input,"Pcdw")
        Pcdw = zeros(Float64,nlayer)
    else
        Pcdw = [parse(Float64,x) for x in split(split(input[indx],"=")[2])]
    end
    #g_tilde
    indx = findfirst(x->contains(x,"g_tilde"),input)
    if indx==nothing
        push!(missing_input,"g_tilde")
        g_tilde = zeros(Float64,nlayer)
    else
        g_tilde = [parse(Float64,x) for x in split(split(input[indx],"=")[2])]
    end
    #t_orth
    indx = findfirst(x->contains(x,"t_orth"),input)
    if indx==nothing
        push!(missing_input,"t_orth")
        t_orth = zeros(Float64,nlayer-1)
    else
        t_orth = [parse(Float64,x) for x in split(split(input[indx],"=")[2])]
    end
    #ac
    indx = findfirst(x->contains(x,"ac"),input)
    if indx==nothing
        push!(missing_input,"ac")
        ac = ones(Float64,nlayer-1)
    else
        ac = [parse(Float64,x) for x in split(split(input[indx],"=")[2])]
    end

    ##If we're missing some layer parameters input
    if length(missing_input)>0
        println("\n###############################################################################")
        println("WARNING : the following layer parameters were not defined :",missing_input)
        println()
        println("These parameters are set to their default values, which are defined below :
J = 0.0\nchi = 0.0\ndoping = 1.0\nD = 0.0\nyrz_mode = Norm\nmu_phys = 0.0#
mu_aux  = 0.0\nQ  = []\nPcdw = 0.0\ng_tilde = 0.0\nt_orth = 0.0\nac = 1.0")
        println("###############################################################################\n")
    end

    pLayer = LayerParameters(nlayer=nlayer,t0=t0,tp0=tp0,Q=Q,Pcdw=Pcdw,D=D,yrz_mode=yrz_mode,mu=mu,mu_aux=mu_aux,x=x,chi=chi,J=J,t_orth=t_orth,ac=ac)


    return pSim,pLayer
end
