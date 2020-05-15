import sys
import time

# Import the MDI Library
try:
    import mdi.MDI_Library as mdi
except:
    raise Exception("Unable to import the MDI Library")

# Import MPI Library
#try:
#    from mpi4py import MPI
#    use_mpi4py = True
#    mpi_comm_world = MPI.COMM_WORLD
#except ImportError:
use_mpi4py = False
mpi_comm_world = None

import math
import utils.collectivevariable as cv
import utils.distance as distance
import utils.utils as ut
try:
    import utils.plot as pl
    animated_draw = True
except ImportError:
    animated_draw = False

if __name__ == "__main__":

    print("AAA HERE 1")

    # Read the command-line options
    iarg = 1
    mdi_options = None
    while iarg < len(sys.argv):
        arg = sys.argv[iarg]

        if arg == "-mdi":
            mdi_options = sys.argv[iarg + 1]
            iarg += 1
        else:
            raise Exception("Unrecognized command-line option")

        iarg += 1

    print("AAA HERE 2")
    # Confirm that the MDI options were provided
    if mdi_options is None:
        raise Exception("-mdi command-line option was not provided")

    # Initialize the MDI Library
    mdi.MDI_Init(mdi_options, mpi_comm_world)

    print("AAA HERE 3")
    # Get unit conversions
    colvar = distance.Distance(319, 320)
    kcalmol_to_atomic = mdi.MDI_Conversion_Factor("kilocalorie_per_mol","atomic_unit_of_energy")
    angstrom_to_atomic = mdi.MDI_Conversion_Factor("angstrom","atomic_unit_of_length")
    kcalmol_per_angstrom_to_atomic = kcalmol_to_atomic / angstrom_to_atomic

    # Input parameters
    width = 0.2 * angstrom_to_atomic # Gaussian width of first collective variable
    height = 0.1 * kcalmol_to_atomic # Gaussian height of first collective variable
    #total_steps = 30000000 # Number of MD iterations. Note timestep = 2fs
    total_steps = 10000 # Number of MD iterations. Note timestep = 2fs
    tau_gaussian = 400 # Frequency of addition of Gaussians
    upper_restraint = 14.0 * angstrom_to_atomic
    lower_restraint = 1.0 * angstrom_to_atomic
    upper_window = 8.0 * angstrom_to_atomic
    lower_window = 2.4 * angstrom_to_atomic
    k_restraint = 10 * kcalmol_per_angstrom_to_atomic
    verbose = False
    #animated_draw = True

    grid_fac = 1
    ngrid = 4000 * grid_fac
    dgrid = 0.005 * angstrom_to_atomic / float(grid_fac)
    bias = [ 0.0 for i in range(ngrid) ]
    bias_derv = [ 0.0 for i in range(ngrid) ]

    # Timings variables
    time_iter = 0.0
    time_integrate = 0.0
    time_bias_pot = 0.0
    time_bias_force = 0.0
    time_draw = 0.0
    time_force_update = 0.0

    s_of_t = [ ] # value of collective variable at time t'
    print("AAA HERE 4")

    # Creat a plot of the results
    if animated_draw:
        my_plot = pl.AnimatedPlot(kcalmol_to_atomic, angstrom_to_atomic)

    # Connect to the engines
    comm = mdi.MDI_Accept_Communicator()

    # Perform the simulation
    mdi.MDI_Send_Command("<NAME", comm)
    name = mdi.MDI_Recv(mdi.MDI_NAME_LENGTH, mdi.MDI_CHAR, comm)
    print("ENGINE NAME: " + str(name))

    # Get the number of atoms
    mdi.MDI_Send_Command("<NATOMS", comm)
    natoms = mdi.MDI_Recv(1, mdi.MDI_INT, comm)

    # Get the number of masses
    mdi.MDI_Send_Command("<MASSES", comm)
    masses = mdi.MDI_Recv(natoms, mdi.MDI_DOUBLE, comm)

    # Initialize an MD simulation
    mdi.MDI_Send_Command("@INIT_MD", comm)

    print("MD Simulation successfully initialized")

    # Open file for Gaussian outputs
    output = open("s_of_t.out", "w")

    for time_step in range(total_steps + 1):
        time_iter_start = time.clock()

        # Proceed to the next point of force evaluation
        time_start = time.clock()
        mdi.MDI_Send_Command("@FORCES", comm)

        # Get simulation box size
        mdi.MDI_Send_Command("<CELL", comm)
        cell_size = mdi.MDI_Recv(9, mdi.MDI_DOUBLE, comm)

        # Note: The following assumes that the cell vectors are orthogonal
        box_len = [ cell_size[ 4 * i ] for i in range(3) ]

        # Get current Cartesian coordinates
        mdi.MDI_Send_Command("<COORDS", comm)
        coords = mdi.MDI_Recv(3*natoms, mdi.MDI_DOUBLE, comm)
        time_end = time.clock()
        time_integrate += time_end - time_start
        
        # Compute values of CVs and gradients
        time_start = time.clock()
        colvar.Evaluate(coords, natoms, box_len)
        colvar_val = colvar.GetValue()

        # Update the bias function
        if time_step % tau_gaussian == 0:
            for i in range(ngrid):
                arg = ( i * dgrid ) - colvar_val
                bias[i] -= ut.Gaussian(arg, width, height)
                bias_derv[i] += ut.Gaussian_derv(arg, width, height)
        time_end = time.clock()
        time_bias_pot += time_end - time_start

        if time_step % tau_gaussian == 0:
            time_start = time.clock()
            if animated_draw:
                my_plot.show(s_of_t, width, height, bias, dgrid)
            time_end = time.clock()
            time_draw += time_end - time_start

        # Evaluate the derivative of Gaussians wrt to Cartesian Coordinates
        time_start = time.clock()
        index1 = int( math.floor( colvar_val / dgrid ) )
        index2 = index1 + 1
        fgrid = ( colvar_val / dgrid ) - float(index1)
        test = (1.0 - fgrid) * bias_derv[ index1 ] + fgrid * bias_derv[ index2 ]
        #print( "Err: " + str(test) + " " + str( dVg_ds ) + " " + str( (test - dVg_ds) / max( abs(dVg_ds), 0.00000000000001 ) ) )
        dVg_ds = test
        ds_dr = colvar.GetGradient()

        # Apply restraints
        if colvar_val > upper_restraint:
            dVg_ds = k_restraint * ( colvar_val - upper_restraint )
        time_end = time.clock()
        time_bias_force += time_end - time_start

        # Compute the updated forces
        time_start = time.clock()
        mdi.MDI_Send_Command("<FORCES", comm)
        forces = mdi.MDI_Recv(3*natoms, mdi.MDI_DOUBLE, comm)
        for idx_atom in range(2):
            for idx_dir in range(3):
                iatom = colvar.GetAtoms()[idx_atom]
                forces[ 3 * iatom + idx_dir ] -= dVg_ds * ds_dr[idx_atom][idx_dir]

        # Send the new forces to the engine
        mdi.MDI_Send_Command(">FORCES", comm)
        mdi.MDI_Send(forces, 3*natoms, mdi.MDI_DOUBLE, comm)
        time_end = time.clock()
        time_force_update += time_end - time_start

        time_iter_end = time.clock()
        time_iter += time_iter_end - time_iter_start

        if time_step % tau_gaussian == 0:
            time_integrate *= 100.0 / time_iter
            time_bias_pot *= 100.0 / time_iter
            time_bias_force *= 100.0 / time_iter
            time_draw *= 100.0 / time_iter
            time_force_update *= 100.0 / time_iter
            print("Iteration " + str(time_step) + " out of " + str(total_steps) + "   (" + str(time_iter) + " s)")
            print("   Integration:  " + str(time_integrate) + "%")
            print("   Bias Pot:     " + str(time_bias_pot) + "%")
            print("   Bias Force:   " + str(time_bias_force) + "%")
            print("   Draw:         " + str(time_draw) + "%")
            print("   Force Update: " + str(time_force_update) + "%")
            time_integrate = 0.0
            time_bias_pot = 0.0
            time_bias_force = 0.0
            time_draw = 0.0
            time_force_update = 0.0

            time_iter = 0.0

            #print("   Colvar_val: " + str(colvar_val/angstrom_to_atomic))

            output.write(str(i) + " " + str(colvar_val) + " " + str(width) + " " + str(height) + "\n")

    # Send the "EXIT" command to each of the engines
    mdi.MDI_Send_Command("EXIT", comm)
    
    if animated_draw:
        my_plot.finalize()
