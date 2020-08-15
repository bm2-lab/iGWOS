uCRISPR - Unified energetics analysis to evaluate the Cas9 on-target activity and off-target effects.
Version 0.1
Author: Dong Zhang, Travis Hurst, Dongsheng Duan & Shi-Jie Chen
Email: chenshi@missouri.edu
Date: Feb 10, 2019

1. System Requirement

        Linux/Unix
        gcc compiler >4.8 version 

2. Compilation
  
    Step 1. Install RNAstructure package (http://rna.urmc.rochester.edu/RNAstructure.html)
            ***Specify the location of the thermodynamic parameters before the next step
            ***(https://rna.urmc.rochester.edu/Text/Thermodynamics.html)

    Step 2. Set path for DATA (line #12) and RNAstructure (line #13) in "uCRISPR.cpp" file:
            string DATA="[Absolute path for directory in which uCRISPR resides]/uCRISPR/data/";
            string RNAstructure="[Absolute path for directory in which RNAstructure resides]/RNAstructure/exe/EnsembleEnergy"

    Step 3. Make: g++ -o uCRISPR uCRISPR.cpp -std=c++11

3. Usage: uCRISPR [options] 

     Options:
       -h         #Display help information
       -on file   #Evaluate on-target avtivities for on-target sites in file,
                   one site (23-mer or 30-mer sequence) per line.
                   See "./example/OntargetSite.dat" for an example.
       -off file  #Evaluate off-target efficiencies for off-target sites in file,
                   one site (20-mer/23-mer wild type sequence and 23-mer off target
                   sequence) per line in file. See "./example/OfftargetSite.dat" for an example.
    
4. Examples 

     a. uCRISPR -on  ./example/OntargetSite.dat   #Evaluate on-target sites given in file "./example/OntargetSite.dat",
                                                   predicted scores are shown on the screen(Also see ./example/OntargetSite.out).

     b. uCRISPR -off ./example/OfftargetSite.dat  #Evaluate off-target sites given in file "./example/OfftargetSite.dat"
                                                   predicted scores are shown on the screen(Also see ./example/OfftargetSite.out).
