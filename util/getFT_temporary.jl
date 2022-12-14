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

alpha=1.0
Bunits = 27800
nB_FT = 5000
Fmax = 5000
FTmesh = LinRange(0,Fmax,nB_FT)
path="/home/bacb3201/Documents/Projects/QO_TB/Phenomenological_model/benchmark_Mei2016/nx10000/Q2pi4_results/"
file=path*"benchmark_Mei_nx100000_eta0.0005.DOS_P000"
w1, dos1 = read_dos(file)
fit1 = fit(w1,dos1,0)                                       #Use polynomial fit of order n for background

path="/home/bacb3201/Documents/Projects/QO_TB/Phenomenological_model/benchmark_Mei2016/nx10000/P003_new/"
file=path*"benchmark_Mei_nx100000_eta0.001.DOS"
w2, dos2 = read_dos(file)
fit2 = fit(w2,dos2,0)                                       #Use polynomial fit of order n for background

for ipts=1:length(w1)
    dos1[ipts] = (dos1[ipts]-fit1(w1[ipts]))*besseli(1,pi*alpha*sqrt(1-((2*(ipts-1))/(length(w1)-1)-1)^2))/besseli(1,pi*alpha)
 #   dos2[ipts] = (dos2[ipts]-fit2(w2[ipts]))*besseli(1,pi*alpha*sqrt(1-((2*(ipts-1))/(length(w2)-1)-1)^2))/besseli(1,pi*alpha)
end

F1 = do_FT(nB_FT,FTmesh,w1./Bunits,dos1)
#F2 = do_FT(nB_FT,FTmesh,w2./Bunits,dos2) #fft(dos)
#freqs = fftfreq(length(w), 1.0/((-w[1]+w[length(w)])/length(w)))

dos_domain = plot((1. ./w1)*Bunits, [dos1], title = "dos")
#dos_domain = plot(fit1,extrema(w1)...,label="fit")
spec_domain = plot(FTmesh,[abs.(F1)],xlim=(0,2000), title = "Spectrum",xticks=[0,200,400,600,800,1000,1200,1400,1600,1800,2000])
spec_domain = vline!([835])
#spec_domain = vline!([640*i for i=1:5])
#spec_domain = vline!([1330*i for i=1:5])
#spec_domain = vline!([440])
#spec_domain = vline!([2*640-170])
#plot(dos_domain, layout = (1,1), grid=true)
plot(dos_domain, spec_domain, layout = (2,1), grid=true)