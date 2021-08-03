all_dat <- read.csv("All_Rat_Data_final_06_01_2021.csv", row.names = 1)
animaln = "A33"

data_t = c(0.3,10,24,48,72,168) #rounded 0.3 hours to 1, removed areas where model converts data_t from days to hours.
data_v = c(all_dat$Volume_Norm[which(all_dat$Animal == animaln & all_dat$Time == 0.3)],
         all_dat$Volume_Norm[which(all_dat$Animal == animaln & all_dat$Time == 10)],
         all_dat$Volume_Norm[which(all_dat$Animal == animaln & all_dat$Time == 24)],
         all_dat$Volume_Norm[which(all_dat$Animal == animaln & all_dat$Time == 48)],
         all_dat$Volume_Norm[which(all_dat$Animal == animaln & all_dat$Time == 72)],
         all_dat$Volume_Norm[which(all_dat$Animal == animaln & all_dat$Time == 168)]) #,
#all_dat$Volume_Norm[which(all_dat$Animal == animaln & all_dat$Time == 336)])
#another matrix for Ultrasound data. PV_flow, Stiffness, and Oxygenation are the rows in order. Columns are time, same as above
# dat_us <- matrix(unlist(c(as.vector(all_dat[which(all_dat$Animal == animaln & all_dat$Time == 0.3), c("PV_flow", "Stiffness", "Oxygenation")]),
#                           as.vector(all_dat[which(all_dat$Animal == animaln & all_dat$Time == 10), c("PV_flow", "Stiffness", "Oxygenation")]),
#                           as.vector(all_dat[which(all_dat$Animal == animaln & all_dat$Time == 24), c("PV_flow", "Stiffness", "Oxygenation")]),
#                           as.vector(all_dat[which(all_dat$Animal == animaln & all_dat$Time == 48), c("PV_flow", "Stiffness", "Oxygenation")]),
#                           as.vector(all_dat[which(all_dat$Animal == animaln & all_dat$Time == 72), c("PV_flow", "Stiffness", "Oxygenation")]),
#                           as.vector(all_dat[which(all_dat$Animal == animaln & all_dat$Time == 168), c("PV_flow", "Stiffness", "Oxygenation")]))),
#                  #as.vector(all_dat[which(all_dat$Animal == animaln & all_dat$Time == 336), c("PV_flow", "Stiffness", "Oxygenation")]))),
#                  ncol = length(data_v))
# 
# dat_stack <- rbind(data_v, dat_us)
# dat_stack[4,] <- dat_stack[4,] / 100

initialFraction = data_v[1] # Fraction of liver mass remaining (N)
#initialFraction = data[1,1]

#Setting initial x-values for ODE solver
cellfractions = c(initialFraction,0,0)
molecularconc = ones(1,7)
G0 = 1

#added US variables to x vector of ODEs
#Pre_dat <- matrix(unlist(c(as.vector(all_dat[which(all_dat$Animal == animaln & all_dat$Time == -24), c("PV_flow", "Stiffness", "Oxygenation")]))))

#Portal_Flow = 0.75 # m/(cm^2*s)
#Stiffness = 2.0 # m/s elastography
#Oxygenation = 0.5 # percent Hb

#x0 = c(cellfractions, molecularconc, G0)
x_0 = as.matrix(c(cellfractions, molecularconc, G0
	#, Portal_Flow, Stiffness, Oxygenation
	))
