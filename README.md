# CPT-singletop

1. under CPT-singletop, mkdir results, mkdir data
2. under results, mkdir (uncertainty, modulation, matrix_elements, graph, histo)
3. 1Mt13TeV.root,1Mtbar13TeV.root (needed only for MG0) can be downloaded here: https://drive.google.com/drive/folders/1Ok0V6b0QI3uDAX2IOOsWYCSkdu8jpbPK?usp=sharing
    MG1.root, MG2.root, MG3.root are on lyoserv. Go to /commun/MG5_aMC_v2_7_2/ExRootAnalysis
4. under data, put 1Mt13TeV.root,1Mtbar13TeV.root, MG1.root, MG2.root, MG3.root 

5. setting.cpp is under CPT-singletop/include
  - Description : It sets some parameters for each case (MG0, MG1, MG2, MG3)
  - Compilation: If you change any things in this program, you have to do "make clean" before "make" 
  - It will ask which simulation you want to analyse when you do step 7, 8, 9, 10

6. Compilation : in CPT-singletop, do "make clean" then "make" 
7. Run b_mu2020.cpp twice, once with t, and once with tbar
8. Run Sensibility.cpp
9. Run XCarre.cpp
10. Run syst_error.cpp
