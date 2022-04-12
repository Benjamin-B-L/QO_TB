include("params.jl")
include("layers.jl")
include("bareFS.jl")
include("dos_wmagfield.jl")
include("IO.jl")


##Get the input file
if length(ARGS)!=1
    throw(ArgumentError("Number of arguments should be equal to one and argument should refer to input file"))
end
inputfile = ARGS[1]

pSim,pLayer,BareFSswitch,DOSswitch = initialize(inputfile)


if BareFSswitch
    #Compute the bare Fermi Surface & write it to file
    FS = getBareFS(pLayer,pSim.nk,0.01)
    saveFS(FS,pSim)
end

if DOSswitch
    print("bla")
    #Compute DOS for all magnetic field values
    dos = getBdepDOS(pSim,pLayer)
    saveDOS(dos,pSim)
end
