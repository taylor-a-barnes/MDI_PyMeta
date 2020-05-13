import sys

# Import the MDI Library
try:
    import mdi.MDI_Library as mdi
except:
    raise Exception("Unable to import the MDI Library")

# Import MPI Library
try:
    from mpi4py import MPI
    use_mpi4py = True
    mpi_comm_world = MPI.COMM_WORLD
except ImportError:
    use_mpi4py = False
    mpi_comm_world = None

import utils.collectivevariable as cv
import utils.distance as distance
import utils.utils as ut
import utils.plot as pl

if __name__ == "__main__":

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

    # Confirm that the MDI options were provided
    if mdi_options is None:
        raise Exception("-mdi command-line option was not provided")

    # Initialize the MDI Library
    mdi.MDI_Init(mdi_options, mpi_comm_world)

    # Get unit conversions
    colvar = distance.Distance(319, 320)
    kcalmol_to_atomic = mdi.MDI_Conversion_Factor("kilocalorie_per_mol","atomic_unit_of_energy")
    angstrom_to_atomic = mdi.MDI_Conversion_Factor("angstrom","atomic_unit_of_length")
    kcalmol_per_angstrom_to_atomic = kcalmol_to_atomic / angstrom_to_atomic

    # Input parameters
    width = 0.2 * angstrom_to_atomic # Gaussian width of first collective variable
    height = 0.1 * kcalmol_to_atomic # Gaussian height of first collective variable
    total_steps = 8000 # Number of MD iterations. Note timestep = 2fs
    tau_gaussian = 400 # Frequency of addition of Gaussians
    upper_restraint = 14.0 * angstrom_to_atomic
    lower_restraint = 1.0 * angstrom_to_atomic
    upper_window = 8.0 * angstrom_to_atomic
    lower_window = 2.4 * angstrom_to_atomic
    k_restraint = 10 * kcalmol_per_angstrom_to_atomic
    verbose = False;

    s_of_t = [ ] # value of collective variable at time t'

    # Creat a plot of the results
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

    for time_step in range(total_steps + 1):
        print("Iteration " + str(time_step) + " out of " + str(total_steps))

        # Proceed to the next point of force evaluation
        mdi.MDI_Send_Command("@FORCES", comm)

        # Get simulation box size
        mdi.MDI_Send_Command("<CELL", comm)
        cell_size = mdi.MDI_Recv(9, mdi.MDI_DOUBLE, comm)

        # Note: The following assumes that the cell vectors are orthogonal
        box_len = [ cell_size[ 4 * i ] for i in range(3) ]

        # Get current Cartesian coordinates
        mdi.MDI_Send_Command("<COORDS", comm)
        coords = mdi.MDI_Recv(3*natoms, mdi.MDI_DOUBLE, comm)

        # Compute values of CVs and gradients
        colvar.Evaluate(coords, natoms, box_len)
        colvar_val = colvar.GetValue()

        # Update the bias function
        if time_step % tau_gaussian == 0:
            s_of_t.append(colvar_val)
            my_plot.show(s_of_t, width, height)

        # Evaluate the derivative of Gaussians wrt to Cartesian Coordinates
        dVg_ds = 0.0
        for gauss in s_of_t:
            s_of_x = colvar_val
            arg = s_of_x - gauss
            dVg_ds += ut.Gaussian_derv(arg, width, height)
        ds_dr = colvar.GetGradient()

        # Apply restraints
        if colvar_val > upper_restraint:
            dVg_ds = k_restraint * ( colvar_val - upper_restraint )

        # Compute the updated forces
        mdi.MDI_Send_Command("<FORCES", comm)
        forces = mdi.MDI_Recv(3*natoms, mdi.MDI_DOUBLE, comm)
        for idx_atom in range(2):
            for idx_dir in range(3):
                iatom = colvar.GetAtoms()[idx_atom]
                forces[ 3 * iatom + idx_dir ] -= dVg_ds * ds_dr[idx_atom][idx_dir]

        # Send the new forces to the engine
        mdi.MDI_Send_Command(">FORCES", comm)
        mdi.MDI_Send(forces, 3*natoms, mdi.MDI_DOUBLE, comm)

        print("   Colvar_val: " + str(colvar_val/angstrom_to_atomic))

    # Send the "EXIT" command to each of the engines
    mdi.MDI_Send_Command("EXIT", comm)

    # Print the data to an output file
    output = open("s_of_t.out", "w")
    for i in range(len(s_of_t)):
        output.write(str(i) + " " + str(s_of_t[i]) + "\n")
    
    my_plot.finalize()
