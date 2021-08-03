#Parameter Sweep on M value
library(knitr) 
library(fOptions)
library(pracma)
library(NlcOptim)
library(openxlsx)
#library(deSolve)
library(sensitivity)
library(parallel)
library(ggplot2)
library(dplyr)
library(stats)
library(parallel)
library(doParallel)

rm(list = ls())

#Set this to the unzipped folder
setwd("~/")


k0_guess = c()
k0_guess[1] = 20.8 # 10.041128 # Metabolic load (M) k(2) Rat - 20.8
k0_guess[2] = 1.5 # TNF production k(3)
k0_guess[3] = 0.9 # TNF mRNA degrad k(4)
k0_guess[4] = 2e4 # JAK activation rate (V_JAK)
k0_guess[5] = 1e4 # Km JAK
k0_guess[6] = 0.4 # JAK degrad
k0_guess[7] = 2 # Conc. of monomeric STAT3
k0_guess[8] = 7.5e2 # STAT3 activation rate
k0_guess[9] = 0.4 # Km STAT3
k0_guess[10] = 0.1 # STAT3 degrad
k0_guess[11] = 2.4e4 # SOCS activation rate
k0_guess[12] = 7e-4 # Km SOCS3
k0_guess[13] = 0.4 # SOCS3 degrad.
k0_guess[14] = 1.5e-2 # SOCS3 inhibition constant
k0_guess[15] = 2.5e2 # IE response gene expression rate
k0_guess[16] = 18 # Km for IE genes
k0_guess[17] = 5 # IE gene degrad
k0_guess[18] = 7 # ECM degrad by MMPs
k0_guess[19] = 33 # ECM degrad

# Cellular parameters
k0_guess[20] = 0.113 # GF production
k0_guess[21] = 0.23 # GF degrad
k0_guess[22] = 6e-2 # k_up - GF uptake/production rate by ECM
k0_guess[23] = 7e-3 # k_Q
k0_guess[24] = 4.4e-3 # k_P
k0_guess[25] = 5.4e-2 # k_R ##Rat - 5.4e-2 # previously (incorrectly) using 5.4e-3
k0_guess[26] = 2e-2 # k_prol - specifies length of mitotic cycle
k0_guess[27] = 0.1 # k_req - requiescence rate
k0_guess[28] = 8 # theta_req (for sigma_req)
k0_guess[29] = 3 # beta_req (for sigma_req)
k0_guess[30] = 0.1 #1e-2 #1e-2 # k_cd - cell death rate ##Rat - 0.1
k0_guess[31] = 9e-3 #0.011050 # theta_cd #Rat 9e-3
k0_guess[32] = 4.5e-3 #0.026046 # beta_cd #Rat 4.5e-3
k0_guess[33] = 3.5e-4 #6.5675e-4# growth constant for cell mass, formerly k(41) #Rat - 3.5e-4

#added ultrasound parameters, test two
#k0_guess[34] = 1    #  portal flow rate increase by Regen Stim
#k0_guess[35] = 1      # Stiffness Increase by PF
#k0_guess[36] = 1     # Stiffness self-reg
#k0_guess[37] = 1    # Oxygenation increase by PF
#k0_guess[38] = 1   # Oxygenation auto-reg

lower <- k0_guess / 3
upper <- k0_guess * 3
lower[30:33] <- k0_guess[30:33] / 1.2
upper[31:33] <- k0_guess[31:33] * 1.2
upper[30] <- k0_guess[30]*1.05


eps = 0.01
Nss = 0.99


#times = seq(from = tStart, to = tEnd, by = 0.1)
t_0 = 0
t_f = 400 #90*24 # 910 * 24 = ~2.5 years


LiverRegenODE <- function(t,x,k){
  #browser()
  # Read-in Tunable Parameters: Disabled to set option in code chunk where LiverRegenODE is called
  #k <- as.matrix[k0_guess]
  
  #Set variable names
  G = x[11];
  Q = max(0,x[1]); P = max(0,x[2]); R = max(0,x[3]); 
  #N_0 = 1
  N = (Q+G*(P+R))+eps; # Calculate N = Q + G(P+R) + epsilon
  Npre = 1+eps;
  #print(Npre)
  IE = x[4]; GF = max(1,x[5]); ECM = x[6]; TNF = x[7]; JAK = x[8]; STAT = x[9]; SOCS = x[10];
  IE0 = x_0[4]; GF0 = x_0[5]; ECM0 = x_0[6]; TNF0 = x_0[7]; JAK0 = x_0[8]; STAT0 = x_0[9]; SOCS0 = x_0[10];
  
  # PF = x[12]; Stiff = x[13]; Oxy = x[14];
  # PF0 = x_0[12]; Stiff0 = x_0[13]; Oxy0 = x_0[14]
  
  # Calculate sigmas
  sr = 0.5*(1+tanh((k[28]-GF)/k[29])); # sigma_req
  sa = 0.5*(1+tanh((k[31]-((N)/k[1]))/k[32])); # sigma_ap
  
  # Steady-State Constants from tunable Parameters
  
  ss1 = -k[2]*k[1] + k[4]/(1+k[5]) + k[3] #-22.3;  TNF production, formerly k34
  ss2 = -k[4]/(1+k[5]) + k[6] #-1.6; JAK production, formerly k35
  ss3 = -k[8]*k[7]^2/(k[7]^2+k[9]*(1+1/k[14])) + 
    k[15]/(1+k[16]) + k[11]/(1+k[12]) + k[10] #2.39e4; % STAT3 production, formerly k36
  ss4 = -k[11]/(1+k[12]) + k[13] #2.4e4; % SOCS production, formerly k37
  ss5 = -k[15]/(1+k[16]) + k[17] #-8.16; % IE gene production, formerly k38
  ss6 = k[18] + k[19] #40; % ECM production, formerly k39
  ss7 = -k[20]*k[1]/(Nss+eps) + k[22] + k[21] #-1.6; % GF production, formerly k40
  
  # Calculate reaction rates
  r1 = k[23]*(IE-IE0)*Q; # Priming of Quiescent cells
  r2 = k[25]*ECM*R; # Replicating cell returning to Quiescence
  r3 = k[27]*sr*P; # Requiescence of Primed cells
  r4 = k[24]*(GF-GF0)*P; # Primed cells begining Replication
  r5 = k[26]*R; # Doubling of Replicating cells
  ra1 = k[30]*sa*Q; # cell death of Quiescent cells
  ra2 = k[30]*sa*P; # cell death of Primed cells
  ra3 = k[30]*sa*R; # cell death of Replicating cells
  
  r6 = k[2]*k[1]/(N); # TNF production by stimulus
  r7 = k[4]*TNF/(TNF+k[5]); # TNF production of JAK
  r8 = k[3]*TNF; # TNF degradation
  r9 = k[6]*JAK; # JAK degredation
  r10 = k[8]*JAK*k[7]^2/(k[7]^2+k[9]*(1+SOCS/k[14])); # STAT3 production by JAK
  r11 = k[15]*STAT/(STAT+k[16]) # IE production by STAT3
  r12 = k[11]*STAT/(STAT+k[12]) # SOCS3 production by STAT3
  r13 = k[10]*STAT # STAT3 degredation
  r14 = k[13]*SOCS # SOCS3 degredation
  r15 = k[17]*IE # IE degredation
  r16 = k[18]*TNF*ECM; # ECM degredation by TNF activated MMPs
  r17 = k[19]*ECM; # ECM degredation
  r18 = k[20]*k[1]/N; # GF production by stimulus
  r19 = k[22]*GF*ECM; # GF uptake/production by ECM
  r20 = k[21]*GF; # GF degredation
  #r21 = k[43]*(PF-PF0) # GF release due to increased Portal_flow
  
  dPh = r1 - r4 - r3 - ra2
  dRh = r4 - r2 + r5 - ra3
  #print(k[35] * ((Npre/N)*PF0 - PF) + k[36] * ((Npre/N)*Stiff0 - Stiff))
  as.matrix(c(-r1 + r2 + r3 - ra1, # Q Phase
              dPh, # P Phase
              dRh, # R Phase
              r11 - r15 + ss5, # IE genes
              r18 - r19 - r20 + ss7, # + r21 GF
              -r16 - r17 + ss6, # ECM
              r6 - r7 - r8 + ss1, # TNF
              r7 - r9 + ss2, # JAK
              r10 - r11 - r12 - r13 + ss3, # STAT3
              r12 - r14 + ss4, # SOCS3
              (k[1]/(N))*k[33] - k[1]*k[33]#, # Cell Mass (M)
              
              # k[34] * ((Npre/N)*PF0 - PF), # Portal Flow
              # k[35] * ((Npre/N)*PF0 - PF) - k[36] * ((Npre/N)*Stiff0 - Stiff), # Stiffness
              # k[37] * ((Npre/N)*PF0 - PF) - k[38] * ((Npre/N)*Oxy0 - Oxy) # Oxygenation
              # 
              # k[35] * (PF-PF0) + k[36] * (ECM-ECM0) - k[37] * (Stiff - Stiff0),# Stiffness
              # k[38] * ((Npre/N)*Oxy0 - Oxy) # Oxygenation
              
              # k[34]*(k[1]/N) - k[35]*PF - k[36], # Portal Flow
              # k[37] * ECM + k[38] * (PF - PF0) - k[39]*Stiff, # Stiffness
              # k[40] - k[41] * (N) * Oxy - k[42] * Oxy # Oxygenation
  ))
}

M_vals <- seq(from = 5, to = 30, length.out = 100)


#Unweighted Parameter Range:
#cd_vals <- seq(from = 0.01, to = 1, length.out = 100)
#Weighted Parameter Range
cd_vals <- exp(seq(log(5e-4), log(0.5), length.out = 100))


Run_Scan <- function(ParamN) {
  #browser() #This is only for troubleshooting
  print(paste('Collecting Volume Data for all Values of M and Kcd'))
  
  #Source random animal file to run simulations:
  relRfile <- paste("11171", "_datagather_2021.R", sep = "")
  source(file = as.character(relRfile))
  
  #Initialize empty results dataframe: 
  CombScanVols <- matrix(data = NA, nrow = 100, ncol = 100)
  
  #Iterate through all parameter combos, and store volume at 168 hr time point
  for (j in 1:length(cd_vals)){
    k0_guess[ParamN] <- cd_vals[j]
    print(paste("Progress:", j-1, "% scans complete", sep = " "))
    for(i in 1:length(M_vals)){
      k0_guess[1] <- M_vals[i]
      out_1 <- ode23s(LiverRegenODE, t_0, t_f, x_0, k = as.matrix(k0_guess), hmax = 0.1)
      xf_1 <- out_1$y
      tf_1 <- out_1$t
      Nf = as.matrix(xf_1[,1]+xf_1[,11]*(xf_1[,2]+xf_1[,3]))

      Vol1 <- Nf[which.min(abs(tf_1 - 168))]
      CombScanVols[i,j] <- Vol1
    }
  }
  rownames(CombScanVols) <- paste("M_", round(M_vals, digits = 2), sep = "")
  colnames(CombScanVols) <- paste("Kcd_", round(cd_vals, digits = 5), sep = "")
  return(CombScanVols)
}

#
res_ults <- Run_Scan(ParamN = 30)
saveRDS(res_ults, file = "Volume_Data_weighted_range.RDS")
print("Done! Check that result saved as .RDS file!")
  


# # nCors <- detectCores() / 3
# # res_ults <- mclapply(Pars_to_scan, Run_Scan, mc.cores = nCors)
# 
# saveRDS(res_ults, file = "LR_3p_ver1_result.RDS")

