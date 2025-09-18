# Decision-Programming-new-formulation

Folder "oil" contains the extended oil wildcatter codes.
Folder "turbine" contains the turbine inspection and maintenance codes.
Folder "water_cvar" contains the risk averse water management codes.
Folder "Results" contains the results used in the paper. 
Additionally, the folder contains the seeds used to create random instances in the computational experiments.
Folder "Printing" contains codes to calculate statistics over multiple random instances. Note that
to run the codes within this folder, you must start julia server from this folder. Alternatively, you 
can edit the file paths according to your preferences.

If you wish to reproduce the results of the paper, use the seeds given in the files. Then all the scripts
should work without modifications.
If you wish to test the codes with new random instances, follow this procedure:

1. Navigate to the folder of the example that you want to solve.
2. Open *_DP.jl file. Comment out the lines where seeds are read from files and uncomment lines where new seeds are created.
   These spots are marked with todo. Similarly, uncomment at the end of the script lines where the seeds are saved to a file.
3. First, you have to run the *_DP.jl file.
   It generates the seed for each instance and saves it to a file
4. After the *_.DP.jl has terminated, you may run the other julia scripts.

NOTE: you have to start julia with 2 threads to run the code. You can use the following command line prompt:
julia --threads=2
   
