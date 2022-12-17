using FFTW
using Plots
using SpecialFunctions
using Polynomials

t0 = 0              # Start time 
fs = 44100          # Sampling rate (Hz)
tmax = 0.1          # End time       

t = t0:1/fs:tmax;   
signal = sin.(2π * 60 .* t)

F = fft(signal)
freqs = fftshift(fftfreq(length(t), fs))

# plots 
time_domain = plot(t, signal, title = "Signal", label='f',legend=:top)
freq_domain = plot(abs.(F), title = "Spectrum", xlim=(-100, +100), xticks=-100:20:100, label="abs.(F)",legend=:top) 


function read_dos(filename)
    f = open(filename) 
    lines=readlines(f)
    
    #Read nbr of ky points
    nky = parse(Int64,split(lines[1],"=")[2])
    
    #Get nbr of pts and read everything
    npts = length(lines)-2
    dos = zeros(Float64,npts)
    w = zeros(Float64,npts)
    for ipts = 1:npts
        w[ipts] = parse(Float64,split(lines[ipts+2]," ")[1])
        dos[ipts] = parse(Float64,split(lines[ipts+2]," ")[nky+1])
    end
    return w, dos
end

function do_FT(nB_FT,FTmesh,x,y)
    F = zeros(Complex{Float64},nB_FT)

    for iFT=1:nB_FT
        for ix=1:length(x) 
            F[iFT] += exp(im*2*pi*FTmesh[iFT]*x[ix])*y[ix]/sqrt(length(x))
        end
    end

    return F
end

alpha=0.0001
Bunits = 27800
nB_FT = 5000
Fmax = 10000
FTmesh = LinRange(0,Fmax,nB_FT)
path="/home/bacb3201/Documents/Projects/QO_TB/Phenomenological_model/benchmark_Mei2016/nx10000/Q2pi4_results/"
file=path*"benchmark_Mei_nx100000_eta0.002.DOS_P002"
w1, dos1 = read_dos(file)
fit1 = fit(w1,dos1,0)                                       #Use polynomial fit of order n for background

#path="/home/bacb3201/Documents/Projects/QO_TB/Phenomenological_model/benchmark_Mei2016/nx10000/Q2pi4_results/"
#file=path*"benchmark_Mei_nx100000_eta0.002.DOS_P000"
#w2, dos2 = read_dos(file)
#fit2 = fit(w2,dos2,0)                                       #Use polynomial fit of order n for background

for ipts=1:length(w1)
    dos1[ipts] = (dos1[ipts]-fit1(w1[ipts]))#*besseli(1,pi*alpha*sqrt(1-((2*(ipts-1))/(length(w1)-1)-1)^2))/besseli(1,pi*alpha)
#    dos2[ipts] = (dos2[ipts]-fit2(w2[ipts]))*besseli(1,pi*alpha*sqrt(1-((2*(ipts-1))/(length(w2)-1)-1)^2))/besseli(1,pi*alpha)
end

F1 = do_FT(nB_FT,FTmesh,w1./Bunits,dos1)
F2 = do_FT(nB_FT,FTmesh,w2./Bunits,dos2) #fft(dos)
#freqs = fftfreq(length(w), 1.0/((-w[1]+w[length(w)])/length(w)))

#DOS PLOT
dos_domain = plot((1. ./w1)*Bunits, [dos1],line=(:solid,2),label=false, title = "dos, P=0.02/t")
plot!(xlabel="B (T)",ylabel="ρ(T)")
plot!(xtickfontsize=12, ytickfontsize=12)

#SPECTRUM
spec_domain = plot(FTmesh,[abs.(F1)],xlim=(0,3000),ylim=(0,12),line=(:solid,2),label=false)
plot!(title = "Spectrum, P=0.02/t",xticks=range(0,2500,step=500))
plot!(xlabel="F (T)",ylabel="FT[ρ(T)]")
plot!(xtickfontsize=12, ytickfontsize=12)
spec_domain = vline!([445*i for i=1:6],line=(:dash,2),label=false)
spec_domain = vline!([660*i for i=1:2],line=(:dash,2),label=false)
spec_domain = vline!([2420*i for i=1:2],line=(:dash,2),label=false)
spec_domain = vline!([645-440, 2*645-440],line=(:dash,2),label=false)
spec_domain = vline!([2420-440],line=(:dash,2),label=false)
#spec_domain = vline!([640*i for i=1:5])
#spec_domain = vline!([1330*i for i=1:5])
#spec_domain = vline!([440])
#spec_domain = vline!([2*640-170])

#SHOW ALL PLOTS
plot(dos_domain, spec_domain, layout = grid(2,1))
savefig(path*"P002_dos+FT.pdf")