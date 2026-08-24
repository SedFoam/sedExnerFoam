### Prepare the mesh
The command **BlockMesh** will mesh a box without the dune geometry. After the base mesh has been generated, a command **createGeometry.py** read the base mesh and move its vertices to create the initial dune geometry which has a conic shape. For this command to work the mesh has to be written in asciii or the script won't be able to read it. For that, the entry **writeFormat** must be set to **ascii** when running **blockMesh**. After the dune has been created, the **writeFormat** can be set to back to **binary** if needed.

### Run the simulation
- A simulation without bed motion can first be run by setting the entry **bedMotion** to **off** in the file **constant/bedloadProperties** 
- Once the hydrodynamics has reached a steady state (we recommend 10 seconds of simulation here), the results of the simulation without morphodynamics can be used as initial conditions of a second simulation with **bedMotion** set to **on**.
- Ensure that you initial condition folder contain the file **finiteArea/rigidBed** as it is nopied in the result folders. Then run the command **sedExnerFoam** again.
