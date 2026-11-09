# Decision-Programming-new-formulation

The folder Final codes contains the codes used in the final version of "An efficient mixed-integer linear programming formulation for solving influence diagrams".

The files oil_dp_* contain the codes related to our extended oil wildcatter instance. Seq identifies codes related to the sequential instance presented in the Supplementary material. The files pig_farm* contain the implementations of the pig farm problem in the Supplementary material. Finally, test_turbine_discrete* contains the turbine inspection and maintenance implementation and water_test_cvar contains the risk averse water resource management problem instances. The suffix Gurobi refers to instances that are solved with default gurobi parameters. The suffix RJT refers to the RJT model, DP refers to the decision programming model in Hankimaa et al. 2025 and suffix new refers to the models presented in the manuscript. The file value_of_discretization_paper.jl implements the example in Supplementary material Section 2.3.

The seeds used in the instances solved in the paper are in folder Final codes/Data. Note that you should modify the file path in the corresponding code file to match your local machine if you want to use the generated seeds.
   
