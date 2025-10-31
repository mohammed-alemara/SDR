# University of Al-Qadisiyah, 
# College of Administration and Economics, 
# Department of Statistics, 
# Sparse SDR with Reciprocal Elastic Net, M.S.c., 
# Code R program for a simulation study to fit SMAVE.rEN, 
# compare it with five methods, 
# SMAVE,  SMAVE.EN, SMAVE.AdEN, SMAVE.SCAD, SMAVE.MCP,
# This R code automatically sets the working directory, checks for required files,
# It automatically downloads required packages if they are missing.
# Automatically stores results in the current working directory,
# Folder: results1, results2, .....
# This code works only in the R Studio environment.
# Developed by Muhammad M. Al-amara & Professor Dr. Taher Reisan,
# Scientific Review and Supervision, Professor Dr. Ali J. Alkenani.  

       ###################
       # Clear workspace #
       ###################

cat("\014")
rm(list = ls())
if(!is.null(dev.list())) dev.off()
graphics.off()

       #################################################
       # calling the required packages to run the code #
       #################################################

packages <- c("writexl", 
              "glmnet", 
              "far", 
              "MASS", 
              "dr", 
              "mvtnorm", 
              "lars", 
              "rstudioapi",
              "ggplot2",
              "gcdnet",
              "ncvreg",
              "reshape2",
              "caret",
              "parallel")

install_and_load <- function(package) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package)
  }
  library(package, character.only = TRUE)
  return(paste(package, "loaded successfully"))
}

results <- lapply(packages, install_and_load)
print(results)


       ##################
       # Function calls #
       ##################


current_file_path <- rstudioapi::getActiveDocumentContext()$path

current_dir <- dirname(current_file_path)
cat("Current working directory is:", current_dir, "\n")

setwd(current_dir)

files_in_directory <- list.files(current_dir)
cat("Files in current directory:\n")
print(files_in_directory)

required_files <- c("MAVE.LASSO.R",
                    "MAVE.EN.R",
                    "MAVE.AdLASSO.R",
                    "MAVE.rEN.R",
                    "REN.R",
                    "MAVE.SCAD.R",
                    "MAVE.MCP.R")

for (file in required_files) {
  file_path <- file.path(current_dir, file)
  if (file.exists(file_path)) {
    tryCatch({
      source(file_path)
      cat(paste("Successfully sourced:", file, "\n"))
    }, error = function(e) {
      cat(paste("Error sourcing:", file, "\nError message:", e$message, "\n"))
    })
  } else {
    cat(paste("File not found:", file, "\n"))
  }
}

       #########################
       # set simulation option #
       #########################

n       <- 50 #case 3
p       <- 24
#simulation1
#beta1   <- c(1,1,1,1,rep(1,20))


#simulation2
beta1   <- c(1,1,1,1,rep(0,20))

#simulation3
#beta1   <- c(3,1.5,0,0,2,rep(0,19))

#simulation4
#beta1   <- c(1,1,1,1,rep(0,20))     #case 1
#beta2   <- c(rep(0,20),1,1,1,1)     #case 1

#beta1    <- c(1,1,0,1,0,1,rep(0,18)) #case 2
#beta2    <- c(rep(0,18),1,0,1,0,1,1) #case 2

#beta1   <- c(rep(1,12),rep(0,12))   #case 3
#beta2   <- c(rep(0,12),rep(1,12))   #case 3

#simulation5
#beta1   <- c(1,rep(0,23))       #case 1
#beta2   <- c(0,1,rep(0,22))     #case 1

Mu      <- rep(0, p)  
k       <- 1
rho     <- 0.9

       #########################################################
       # where Σ has one of the following covariance structure #
       #########################################################


case    <- "AR"

# Case I (IS): Isotropic design
if (case == "IS") {
  Sigma <- diag(p)
} else if (case == "CS") {
  # Case II (CS): Compound symmetry design
  Sigma <- matrix(rho, nrow = p, ncol = p)
  diag(Sigma) <- 1
} else if (case == "AR") {
  # Case III (AR): Autoregressive design
  Sigma <- outer(1:p, 1:p, function(i, j) rho^abs(i - j))
}

       ################################################
       # The models used in the simulation experiment#
       ###############################################
model.selection <- "model1"

       #####################
       # no. of iterations #
       #####################

iterations     <- 5
start.time     <- Sys.time()

       ####################################
       # Setting up variables for storage #
       ####################################

MSE.MAVE.LASSO   <- numeric(iterations)
MSE.MAVE.EN      <- numeric(iterations)
MSE.MAVE.rEN     <- numeric(iterations)
MSE.MAVE.AdLASSO <- numeric(iterations)
MSE.MAVE.SCAD    <- numeric(iterations)
MSE.MAVE.MCP     <- numeric(iterations)
#######################################
nz.MAVE.LASSO    <- numeric(iterations)
nz.MAVE.EN       <- numeric(iterations)
nz.MAVE.rEN      <- numeric(iterations)
nz.MAVE.AdLASSO  <- numeric(iterations)
nz.MAVE.SCAD     <- numeric(iterations)
nz.MAVE.MCP      <- numeric(iterations)
#######################################
COR.MAVE.LASSO   <- numeric(iterations)
COR.MAVE.EN      <- numeric(iterations)
COR.MAVE.rEN     <- numeric(iterations)
COR.MAVE.AdLASSO <- numeric(iterations)
COR.MAVE.SCAD    <- numeric(iterations)
COR.MAVE.MCP     <- numeric(iterations)

       ###########################################
       # START ITERATUINS  with Serial execution #
       ###########################################

for (i in 1:iterations) {

          xx     <- rmvnorm(n= n , mean = Mu , sigma = Sigma)
          error  <- rnorm(n)
        
        
        if (model.selection == "model1"){
          d = 1
          y <- beta1 %*% t(xx) + 0.5 * error
        } else if (model.selection == "model2"){
          d = 2
          y <- beta1 %*% t(xx) / (0.5 +((beta2 %*% t(xx) + 0.5)^2)) + 0.2 * error
        } else if (model.selection == "model3"){
          d = 2
          y <- sign(beta1 %*% t(xx)) * log(abs(beta2 %*% t(xx)+5)) + 0.2 * error
        }
        
        

        
       #############################
       # fit model with MAVE-LASSO #
       #############################

bhat.MAVE.LASSO     <- MAVE.LASSO(data.frame(xx), y, d)

MSE.MAVE.LASSO[i]   <- bhat.MAVE.LASSO[[3]]
nz.MAVE.LASSO[i]    <- bhat.MAVE.LASSO[[5]]
COR.MAVE.LASSO[i]   <- bhat.MAVE.LASSO[[4]]


       ##################################
       # fit model with MAVE-Elatic Net #
       ##################################

bhat.MAVE.EN        <- MAVE.EN(data.frame(xx), y, d)

MSE.MAVE.EN[i]      <- bhat.MAVE.EN[[3]]
nz.MAVE.EN[i]       <- bhat.MAVE.EN[[5]]
COR.MAVE.EN[i]      <- bhat.MAVE.EN[[4]]

       ######################################
       # fit model with MAVE-adiptive lasso #
       ######################################

bhat.MAVE.AdLASSO     <- MAVE.AdLASSO(data.frame(xx), y, d)

MSE.MAVE.AdLASSO[i]   <- bhat.MAVE.AdLASSO[[3]]
nz.MAVE.AdLASSO[i]    <- bhat.MAVE.AdLASSO[[5]]
COR.MAVE.AdLASSO[i]   <- bhat.MAVE.AdLASSO[[4]]

       ############################
       # fit model with MAVE-SCAD #
       ############################

bhat.MAVE.SCAD       <- MAVE.SCAD(data.frame(xx), y, d)

MSE.MAVE.SCAD[i]     <- bhat.MAVE.SCAD[[3]]
nz.MAVE.SCAD[i]      <- bhat.MAVE.SCAD[[5]]
COR.MAVE.SCAD[i]     <- bhat.MAVE.SCAD[[4]]

       ###########################
       # fit model with MAVE-MCP #
       ###########################
bhat.MAVE.MCP       <- MAVE.MCP(data.frame(xx), y, d)

MSE.MAVE.MCP[i]     <- bhat.MAVE.MCP[[3]]
nz.MAVE.MCP[i]      <- bhat.MAVE.MCP[[5]]
COR.MAVE.MCP[i]     <- bhat.MAVE.MCP[[4]]

       #############################################
       # fit model with MAVE-reciprocal Elatic Net #
       #############################################

bhat.MAVE.rEN       <- MAVE.EN(data.frame(xx), y, d)

MSE.MAVE.rEN[i]     <- bhat.MAVE.rEN[[3]]
nz.MAVE.rEN[i]      <- bhat.MAVE.rEN[[5]]
COR.MAVE.rEN[i]     <- bhat.MAVE.rEN[[4]]


       ################
       # Progress Bar #
       ################

pb = txtProgressBar(min = 0, max = iterations, style = 3)
setTxtProgressBar(pb, i)
}
end.time = Sys.time()

       ##################
       # END ITERATIONS #
       ##################

       ##################################################
       # Calculate iteration time. Format time as HH:MM #
       ##################################################

time.taken        <- end.time - start.time
time.minutes      <- as.integer(ceiling(as.numeric(time.taken, units = "mins")))
time.hours        <- as.integer(time.minutes / 60)
remaining.minutes <- time.minutes %% 60
formatted.time    <- sprintf("%02d:%02d", time.hours, remaining.minutes)
formatted.time    <- as.character(formatted.time)
formatted.time.with.label <- paste("Time", formatted.time)

       ###############################################
       # combine all MSE results to find Average MSE #
       ###############################################

AMSE.MAVE.LASSO   <- mean(MSE.MAVE.LASSO)
AMSE.MAVE.EN      <- mean(MSE.MAVE.EN)
AMSE.MAVE.AdLASSO <- mean(MSE.MAVE.AdLASSO)
AMSE.MAVE.SCAD    <- mean(MSE.MAVE.SCAD)
AMSE.MAVE.MCP     <- mean(MSE.MAVE.MCP)
AMSE.MAVE.rEN     <- mean(MSE.MAVE.rEN)

AMSE <- data.frame(
                   AMSE.MAVE.LASSO   = as.numeric(AMSE.MAVE.LASSO),
                   AMSE.MAVE.EN      = as.numeric(AMSE.MAVE.EN),
                   AMSE.MAVE.AdLASSO = as.numeric(AMSE.MAVE.AdLASSO),
                   AMSE.MAVE.SCAD    = as.numeric(AMSE.MAVE.SCAD),
                   AMSE.MAVE.MCP     = as.numeric(AMSE.MAVE.MCP),
                   AMSE.MAVE.rEN     = as.numeric(AMSE.MAVE.rEN)
                   )

       ##########################################
       # combine all MSE results to find SD MSE #
       ##########################################

SDMSE.MAVE.LASSO   <- sd(MSE.MAVE.LASSO)
SDMSE.MAVE.EN      <- sd(MSE.MAVE.EN)
SDMSE.MAVE.AdLASSO <- sd(MSE.MAVE.AdLASSO)
SDMSE.MAVE.SCAD    <- sd(MSE.MAVE.SCAD)
SDMSE.MAVE.MCP     <- sd(MSE.MAVE.MCP)
SDMSE.MAVE.rEN     <- sd(MSE.MAVE.rEN)

SDMSE <- data.frame(
                   SDMSE.MAVE.LASSO   = as.numeric(SDMSE.MAVE.LASSO),
                   SDMSE.MAVE.EN      = as.numeric(SDMSE.MAVE.EN),
                   SDMSE.MAVE.AdLASSO = as.numeric(SDMSE.MAVE.AdLASSO),
                   SDMSE.MAVE.SCAD    = as.numeric(SDMSE.MAVE.SCAD),
                   SDMSE.MAVE.MCP     = as.numeric(SDMSE.MAVE.MCP),
                   SDMSE.MAVE.rEN     = as.numeric(SDMSE.MAVE.rEN)
                   )

       ############################
       # combine all Av0s results # 
       ############################

Av0s.MAVE.LASSO   <- mean(nz.MAVE.LASSO)
Av0s.MAVE.EN      <- mean(nz.MAVE.EN)
Av0s.MAVE.AdLASSO <- mean(nz.MAVE.AdLASSO)
Av0s.MAVE.SCAD    <- mean(nz.MAVE.SCAD)
Av0s.MAVE.MCP     <- mean(nz.MAVE.MCP)
Av0s.MAVE.rEN     <- mean(nz.MAVE.rEN)
Av0s <- data.frame(
                   Av0s.MAVE.LASSO   = as.numeric(Av0s.MAVE.LASSO),
                   Av0s.MAVE.EN      = as.numeric(Av0s.MAVE.EN),
                   Av0s.MAVE.AdLASSO = as.numeric(Av0s.MAVE.AdLASSO),
                   Av0s.MAVE.SCAD    = as.numeric(Av0s.MAVE.SCAD),
                   Av0s.MAVE.MCP     = as.numeric(Av0s.MAVE.MCP),
                   Av0s.MAVE.rEN     = as.numeric(Av0s.MAVE.rEN)
                   )

       ###############################
       # combine all corl MU results # 
       ###############################

CORMU.MAVE.LASSO   <- mean(COR.MAVE.LASSO)
CORMU.MAVE.EN      <- mean(COR.MAVE.EN)
CORMU.MAVE.AdLASSO <- mean(COR.MAVE.AdLASSO)
CORMU.MAVE.SCAD    <- mean(COR.MAVE.SCAD)
CORMU.MAVE.MCP     <- mean(COR.MAVE.MCP)
CORMU.MAVE.rEN     <- mean(COR.MAVE.rEN)
CORMU <- data.frame(
                    CORMU.MAVE.LASSO   = as.numeric(CORMU.MAVE.LASSO),
                    CORMU.MAVE.EN      = as.numeric(CORMU.MAVE.EN),
                    CORMU.MAVE.AdLASSO = as.numeric(CORMU.MAVE.AdLASSO),
                    CORMU.MAVE.SCAD    = as.numeric(CORMU.MAVE.SCAD),
                    CORMU.MAVE.MCP     = as.numeric(CORMU.MAVE.MCP),
                    CORMU.MAVE.rEN     = as.numeric(CORMU.MAVE.rEN)
                   )  

       ###############################
       # combine all corl SD results # 
       ###############################

CORSD.MAVE.LASSO    <- sd(COR.MAVE.LASSO)
CORSD.MAVE.EN       <- sd(COR.MAVE.EN)
CORSD.MAVE.AdLASSO  <- sd(COR.MAVE.AdLASSO)
CORSD.MAVE.SCAD     <- sd(COR.MAVE.SCAD)
CORSD.MAVE.MCP      <- sd(COR.MAVE.MCP)
CORSD.MAVE.rEN      <- sd(COR.MAVE.rEN)
CORSD <- data.frame(
                    CORSD.MAVE.LASSO   = as.numeric(CORSD.MAVE.LASSO),
                    CORSD.MAVE.EN      = as.numeric(CORSD.MAVE.EN),
                    CORSD.MAVE.AdLASSO = as.numeric(CORSD.MAVE.AdLASSO),
                    CORSD.MAVE.SCAD    = as.numeric(CORSD.MAVE.SCAD),
                    CORSD.MAVE.MCP     = as.numeric(CORSD.MAVE.MCP),
                    CORSD.MAVE.rEN     = as.numeric(CORSD.MAVE.rEN)

)


       #######################
       # Export  all results #
       ####################### 

AMSE_long  <- melt(AMSE)
SDMSE_long <- melt(SDMSE)
Av0s_long  <- melt(Av0s)
CORMU_long <- melt(CORMU)
CORSD_long <- melt(CORSD)

AMSE_long$value  <- round(as.numeric(AMSE_long$value) , 6)
SDMSE_long$value <- round(as.numeric(SDMSE_long$value), 6)
Av0s_long$value  <- round(as.numeric(Av0s_long$value) , 2)
CORMU_long$value <- round(as.numeric(CORMU_long$value), 6)
CORSD_long$value <- round(as.numeric(CORSD_long$value), 6)

AMSE_values  <- AMSE_long$value
SDMSE_values <- SDMSE_long$value
Av0s_values  <- Av0s_long$value
CORMU_values <- CORMU_long$value
CORSD_values <- CORSD_long$value

total <- data.frame(
  AMSE  = AMSE_values,
  SDMSE = SDMSE_values,
  Av0s  = Av0s_values,
  CORMU = CORMU_values,
  CORSD = CORSD_values
)


# Function to get the next available directory name with an incremented number
get_next_directory <- function(base_path, base_name) {
  i <- 1
  repeat {
    dir_name <- paste0(base_name, i)
    dir_path <- file.path(base_path, dir_name)
    if (!dir.exists(dir_path)) {
      return(dir_path)
    }
    i <- i + 1
  }
}

current_dir <- getwd()
new_dir     <- get_next_directory(current_dir, "results")
dir.create(new_dir, recursive = TRUE)

amse_file_path   <- file.path(new_dir, "AMSE.xlsx")
sdmse_file_path  <- file.path(new_dir, "SDMSE.xlsx")
av0s_file_path   <- file.path(new_dir, "Av0s.xlsx")
corMU_file_path  <- file.path(new_dir, "CORMU.xlsx")
corSD_file_path  <- file.path(new_dir, "CORSD.xlsx")
total_file_path  <- file.path(new_dir, "total.xlsx")

save.image(file = file.path(new_dir, "simulation.RData"))
write_xlsx(AMSE_long , path = amse_file_path)
write_xlsx(SDMSE_long, path = sdmse_file_path)
write_xlsx(Av0s_long , path = av0s_file_path)
write_xlsx(CORMU_long, path = corMU_file_path)
write_xlsx(CORSD_long, path = corSD_file_path)
write_xlsx(total     , path = total_file_path)




       ############
       # END CODE #
       ############


formatted.time.with.label
AMSE
SDMSE
Av0s
CORMU
CORSD
