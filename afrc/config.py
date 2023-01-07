# amino acid residues 
AA_list = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']

# R0 parameters for RMS distances as calibrated against FRC simulations
RIJ_RMS_R0 = {'A': 6.5463, 
              'C': 6.2676, 
              'D': 6.3994, 
              'E': 6.2649, 
              'F': 6.2519, 
              'G': 6.1045, 
              'H': 6.2156, 
              'I': 6.4353, 
              'K': 6.3060, 
              'L': 6.2636, 
              'M': 6.3813, 
              'N': 6.2652, 
              'P': 6.4323, 
              'Q': 6.2547, 
              'R': 6.2790, 
              'S': 6.3161, 
              'T': 6.1995, 
              'V': 6.3204, 
              'W': 6.3000, 
              'Y': 6.3188}

# R0 parameters for mean distances as calibrated against FRC simulations
RIJ_R0 = {'A': 6.0381,
          'C': 5.7826,
          'D': 5.911,
          'E': 5.768,
          'F': 5.7612,
          'G': 5.6324,
          'H': 5.7262,
          'I': 5.9361,
          'K': 5.8272,
          'L': 5.7801,
          'M': 5.8894,
          'N': 5.773, 
          'P': 5.9388,
          'Q': 5.7719,
          'R': 5.7921,
          'S': 5.8364,
          'T': 5.7242,
          'V': 5.8409,
          'W': 5.814,
          'Y': 5.8266}


## Note small edits to G and P from homopolymer based on a greater range of callibration curves 
## original unfixed val
## small uniform increase of 0.0005 to correct small deviation at long chains identified
RG_X0 = { 'A': 0.5355,
          'C': 0.5585, 
          'D': 0.5517,
          'E': 0.5563,
          'F': 0.5521,
          'G': 0.5861,
          'H': 0.5595,
          'I': 0.5433,
          'K': 0.5483,
          'L': 0.5555,
          'M': 0.5451,
          'N': 0.5548,
          'P': 0.5549,
          'Q': 0.5567,
          'R': 0.5481,
          'S': 0.5503,
          'T': 0.5645,
          'V': 0.5521,
          'W': 0.5489,
          'Y': 0.5493}



## R0 values for Rg are back calculated by computing Rg vs N curves for homopolymers where Rg is determined as the mean
## value using the parameters identified from the Lhuilier fits (i.e. those in RG_X0) and then fitting for the intercept
## on a log-log plot. In fact, RG_R0 = 1.3858*(1/RG_X0).
RG_R0 = {'A': 2.5866,
         'C': 2.4812,
         'D': 2.5115,
         'E': 2.4909,
         'F': 2.5097,
         'G': 2.3654,
         'H': 2.4768,
         'I': 2.5499,
         'K': 2.5269,
         'L': 2.4945,
         'M': 2.5416,
         'N': 2.4976,
         'P': 2.4972,
         'Q': 2.4892,
         'R': 2.5278,
         'S': 2.5178,
         'T': 2.4551,
         'V': 2.5097,
         'W': 2.5242,
         'Y': 2.5224}


# resolution in Angstroms
P_OF_R_RESOLUTION=0.05

RE_RG_CORRECTION_FACTORS = [0.9490139479, -0.0015426547]

