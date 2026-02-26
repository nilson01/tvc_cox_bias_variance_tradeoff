############################
# Data analysis for TVC Cox paper
# 
# Selukar
# 2026-02-24
############################

source("anaFunctions.R")


# ============================================================
# 6) Load JASA and construct (futime, status, T_treat, D_treat, W)
# ============================================================
data(jasa, package = "survival")

jasa_dat <- jasa %>%
  transmute(
    futime = as.numeric(futime),
    status = as.integer(fustat == 1),
    
    transplant = as.integer(transplant),
    wait_time  = as.numeric(wait.time),
    
    # ---- baseline W (choose what we want W to be) ----
    # Here: surgery as baseline W (0/1-ish)
    W = as.numeric(surgery)
  )


# Treatment switch time:
# If transplanted: switch at wait.time
# Else: set Inf (so Z(t)=0, D(t)=0 always)
jasa_dat$T_treat <- ifelse(jasa_dat$transplant == 1,
                           pmax(jasa_dat$wait_time, 1e-12),
                           Inf)

# Convention: D_treat is fed as x to tt(); tt() makes D(t)=max(0, t-x)
jasa_dat$D_treat <- jasa_dat$T_treat

# Sanity checks
cat("n =", nrow(jasa_dat), "\n")
cat("events =", sum(jasa_dat$status == 1), "\n")
cat("transplanted =", sum(jasa_dat$transplant == 1), "\n")
cat("finite T_treat =", sum(is.finite(jasa_dat$T_treat)), "\n")

# ============================================================
# 7) Run: choose which models we want
#    Options: "Z", "ZD", "ZW", "ZDW"
# ============================================================

chosen_models <- c("Z", "ZD")
optHorizons <- 365.25*c(0.25,0.5,0.75,1,1.5,2,2.5,3)

out <- run_models_and_clsr(
  dat      = jasa_dat,
  models   = chosen_models,
  horizons = optHorizons,        # or provide a numeric vector of horizons
  print_fits = TRUE
)

out$results$Z$fit
out$results$ZD$fit

sum(jasa)

clsrOut <- rbind(t(out$clsr_table[1:8,c("horizon_time","CLSR")]),
                 t(out$clsr_table[9:16,c("CLSR")]))
clsrOut[1,] <- clsrOut[1,]/365.25
rownames(clsrOut) <- c("Years since acceptance","CT","TH")



# =============================================================
# 8) repeat with bootstrap
# ============================================================

set.seed(2026)
numBoots <- 200

bootCT <- bootTH <- matrix(nrow=numBoots,ncol=length(optHorizons))

for (i in 1:numBoots){
  tmpDat <- jasa_dat[sample(nrow(jasa_dat),replace=TRUE),]
  
  tmpOut <- run_models_and_clsr(
    dat      = tmpDat,
    models   = chosen_models,
    horizons = optHorizons,        # or provide a numeric vector of horizons
    print_fits = FALSE
  )
  
  bootCT[i,] <- t(tmpOut$clsr_table[1:8,c("CLSR")])
  bootTH[i,] <- t(tmpOut$clsr_table[9:16,c("CLSR")])
}

clsrLB <- rbind(apply(bootCT,2,quantile,probs=0.025),
                apply(bootTH,2,quantile,probs=0.025)
                )
clsrUB <- rbind(apply(bootCT,2,quantile,probs=0.975),
                apply(bootTH,2,quantile,probs=0.975)
)

round(clsrLB,3)
round(clsrUB,3)

bootOut <- matrix(paste0(round(clsrOut[2:3,],3),
      " (",
      round(clsrLB,3),
      "-",
      round(clsrUB,3),
      ")"),
      nrow=2
)
rownames(bootOut) <- c("CT","TH")
colnames(bootOut) <- optHorizons/365.25

write.csv(bootOut,"tabBoot.csv")
