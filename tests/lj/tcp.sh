#location of required codes
DRIVER_LOC=$(cat ../locations/MDI_PyMeta)
LAMMPS_LOC=$(cat ../locations/LAMMPS)

#remove old files
if [ -d work ]; then 
  rm -r work
fi

#create work directory
cp -r data work
cd work

#launch LAMMPS
${LAMMPS_LOC} -mdi "-role ENGINE -name LAMMPS -method TCP -port 8021 -hostname localhost" -in in.lj > lammps.out &

#launch driver
python ${DRIVER_LOC} -mdi "-role DRIVER -name driver -method TCP -port 8021" &

wait
