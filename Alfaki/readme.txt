################### STANDARD POOLING PROBLEM INSTANCES #########################
These pooling problem instances are used in the experiments section of the paper 
by Alfaki and Haugland (2010). All the instances and the formulations have been 
written in GAMS.

HOW TO RUN THE TEST INSTANCES:
    I. SYSTEM REQUIREMENT:
        1. Windows/Unix.
        2. GAMS modeling system
        3. Global NLP solver (e.g., BARON).
   II. RUN AN INSTANCE:
        1. The default NLP solver is BARON, you can change it in one of the 
           formulation files; pqmodel.gms, tpmodel.gms and stpmodel.gms.
        2. Write the name of the formulation that you want to use in the first
           line of of the file 'xmodel.gms' e.g., if the formulation you want 
           use is STP-formulation, the first line of this file will look like:

           $include stpmodel.gms

        3. Suppose you want solve instance Adhya1, so in the terminal (windows 
           command prompt) write:
           [pooling instance DIR]$ gams Adhya1.gms 

Alfaki, M., and Haugland, D. (2010). Strong formulation for the pooling problem. 
Journal of Global Optimization. doi: 10.1007/s10898-012-9875-6.
################################################################################