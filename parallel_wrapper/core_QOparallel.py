import numpy as np
import multiprocessing as mp
import subprocess


class SimulationParameters:
    """Simple class to combine all parameters relevant for the simulation"""

    def __init__(self, nx, ny_avg, bmin, bmax, nb, theta_in, eta, nk, outfile):
        self.nx = nx
        self.ny_avg = ny_avg
        self.bmin = bmin
        self.bmax = bmax
        self.nb = nb
        self.theta_in = theta_in
        self.eta = eta
        self.nk = nk
        self.outfile = outfile


class LayerParameters:
    """Simple class to combine all paremeters
    relevant for the different layers"""

    def __init__(self, nlayer, t0, tp0, J, chi, x, D, yrz_mode, mu,
                 mu_aux, Q, Pcdw, g_tilde, t_orth, ac):
        self.nlayer = nlayer
        self.t0 = t0
        self.tp0 = tp0
        self.J = J
        self.chi = chi
        self.x = x
        self.D = D
        self.yrz_mode = yrz_mode
        self.mu = mu
        self.mu_aux = mu_aux
        self.Q = Q
        self.Pcdw = Pcdw
        self.g_tilde = g_tilde
        self.t_orth = t_orth
        self.ac = ac


def write_input(file, DOSswitch, BareFSswitch, pSim, pLayer):
    """Write input file for QO DOS calculation"""
    f = open(file, "w")
    # DOS and FS switch
    f.write("##TYPE OF CALCULATION \n")
    if BareFSswitch:
        f.write("BareFSswitch = true\n")
    else:
        f.write("BareFSswitch = false\n")
    if DOSswitch:
        f.write("DOSswitch = true\n")
    else:
        f.write("DOSswitch = false\n")
    f.write("\n")
    # Simulation parameters
    f.write("##SIMULATION PARAMETERS \n")
    f.write("nx = "+str(pSim.nx)+"\n")
    f.write("ny_avg = "+str(pSim.ny_avg)+"\n")
    f.write("bmin = "+str(pSim.bmin)+"\n")
    f.write("bmax = "+str(pSim.bmax)+"\n")
    f.write("nb = "+str(pSim.nb)+"\n")
    f.write("theta_in = "+str(pSim.theta_in)+"\n")
    f.write("broadening = "+str(pSim.eta)+"\n")
    f.write("nk = "+str(pSim.nk)+"\n")
    f.write("outfile = "+str(pSim.outfile)+"\n")
    # Layer parameters
    f.write("##LAYER PARAMETERS \n")
    f.write("nlayer = "+str(pLayer.nlayer)+"\n")
    f.write("t0 = "+" ".join([str(x) for x in pLayer.t0])+"\n")
    f.write("tp0 = "+" ".join([str(x) for x in pLayer.tp0])+"\n")
    f.write("J = "+" ".join([str(x) for x in pLayer.J])+"\n")
    f.write("chi = "+" ".join([str(x) for x in pLayer.chi])+"\n")
    f.write("doping = "+" ".join([str(x) for x in pLayer.x])+"\n")
    f.write("D_yrz = "+" ".join([str(x) for x in pLayer.D])+"\n")
    f.write("yrz_mode = "+" ".join([x for x in pLayer.yrz_mode])+"\n")
    f.write("mu_phys = "+" ".join([str(x) for x in pLayer.mu])+"\n")
    f.write("mu_aux = "+" ".join([str(x) for x in pLayer.mu_aux])+"\n")
    f.write("Q = "+" ".join(["["+",".join(["("+str(x[0])
            + "," + str(x[1])+")"for x in pLayer.Q[q]])+"]" for q in pLayer.Q])+"\n")
    f.write("Pcdw = "+" ".join([str(x) for x in pLayer.Pcdw])+"\n")
    f.write("g_tilde = "+" ".join([str(x) for x in pLayer.g_tilde])+"\n")
    f.write("t_orth = "+" ".join([str(x) for x in pLayer.t_orth])+"\n")
    f.write("ac = "+" ".join([str(x) for x in pLayer.ac])+"\n")

    f.close()

    return True


def combineDOSfiles(file, bgrid, nb, ncol, nfiles):
    # First take header
    header = (open(file+"_0.DOS", "r").readlines())[0:2]
    # Exctract data
    dos = np.zeros([nb, ncol])
    cnt = 0
    for ifile in range(nfiles):
        f = open(file+"_"+str(ifile)+".DOS", "r")
        dat = f.readlines()
        for iline in range(2, len(dat)):
            dos[cnt, :] = [dat[iline].strip().split()[icol+1]
                           for icol in range(ncol)]
            cnt += 1
        f.close()
    # Save data
    f = open(file+".DOS", "w")
    f.write(header[0])
    f.write(header[1])
    for ib in range(nb):
        f.write(str(1/bgrid[ib])+" ")
        for icol in range(ncol):
            f.write(str(dos[ib, icol])+" ")
        f.write("\n")
    f.close()
    # Remove tmp files
    subprocess.run("rm "+file+"_*.DOS", shell=True)


def call_QO(id, bgrid, DOSswitch, BareFSswitch, pSim, pLayer):
    # Define input file
    inputfile = "tmp.input."+str(id)
    # Set the right bgrid for this worker
    pSim.nb = len(bgrid)
    pSim.bmin = bgrid[0]
    pSim.bmax = bgrid[pSim.nb-1]
    pSim.outfile = pSim.outfile+"_"+str(id)
    # Write the input file
    write_input(inputfile, DOSswitch, BareFSswitch, pSim, pLayer)
    # Call the main routines
    subprocess.run("julia ../src/main_QO.jl "+inputfile, shell=True)
    # Clean up input file
    subprocess.run("rm "+inputfile, shell=True)


def call_parallelQO(DOSswitch, BareFSswitch, pSim, pLayer, ncores):
    """main routine managing the parallel call for QO calculation"""
    # Cut magnetic field grid into ncores part
    b_chunks = np.array_split(np.linspace(
        pSim.bmin, pSim.bmax, pSim.nb), ncores)
    # Launch parallel processes
    pool = mp.Pool(processes=ncores)
    tmp = [pool.apply_async(call_QO, args=(
        i, b_chunks[i], DOSswitch, BareFSswitch, pSim, pLayer)) for i in range(ncores)]
    for x in tmp:
        x.get()
    pool.close()
    # Combine all output files & clean up
    # Bare FS
    if BareFSswitch:
        # just take one and remove the others
        subprocess.run("cp "+pSim.outfile+"_0.FS " + pSim.outfile
                       + ".FS && rm "+pSim.outfile+"_*.FS", shell=True)
    # DOS
    if DOSswitch:
        combineDOSfiles(pSim.outfile, np.linspace(
            pSim.bmin, pSim.bmax, pSim.nb), pSim.nb, pSim.ny_avg+1, ncores)
