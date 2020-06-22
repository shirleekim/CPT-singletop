# CPT-singletop

1. under CPT-singletop, mkdir results, mkdir data
2. under results, mkdir (uncertainty, modulation, matrix_element, graph, histo)
3. under data, put 1Mt13TeV.root,1Mtbar13TeV.root, 1Mt13TeV_2020.root, MG1.root, MG2.root, MG3.root 

4. setting.cpp is under CPT-singletop/include
  - Description : It sets some parameters for each case (MG0, MG1, MG2, MG3)
  - Compilation: If you change any things in this program, you have to do "make clean" before "make" 
  - It will ask which simulation you want to analyse 

5. Compilation : in CPT-singletop, do "make clean" then "make" 
6. Run b_mu2020.cpp twice, once with t, and once with tbar
7. Run Sensibility.cpp
8. Run XCarre.cpp
9. Run syst_error.cpp
