
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

    # Connect to the engines

    # Perform the simulation

    # Send the "EXIT" command to each of the engines
