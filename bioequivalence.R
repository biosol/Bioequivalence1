setwd("../Desktop/Bioequivalencia/")
setwd("../../Documents/R/")
### 1. READ RAW DATA FILE ###
data = read.csv("raw_data.csv", sep = ";", stringsAsFactors = FALSE)

# Organise data by subject and treatment
sub = as.vector(data$Subject.No)
treat = as.vector(data$Treatment)
seq = data$Sequence..1.TR..2.RT.

d = data[order(sub,treat),]

# Set columns of d as vectors
subj = as.vector(d$Subject.No)
treatm = as.vector(d$Treatment)
con = as.vector(d$conc)
time = as.numeric(d$Calculated_Real.PK.time..h.)
theo = as.numeric(d$Theorical.PK.time..h.)
treat = as.vector(d$Treatment)
seq = as.vector(d$Sequence..1.TR..2.RT.)
table = data.frame(subj, seq, treatm, time, theo, con)

## LOG-TRANSFORMATION OF CONCENTRATION DATA
clog = log(table$con)
table$conlog = clog

### 2. GRAPHICAL REPRESENTATION OF DATA ###

library(data.table)
table = data.table(table)

sp = split(table, table$sub) #list of data frames for each subject
#Each data frame contains data about the reference and the test for the same indv

## a. REPRESENTATION OF NORMAL AND LOG-TRANSFORMED CONCENTRATIONS VS TIME
pdf("time_conc_plots.pdf")

for (i in sp){
  t = as.vector(i$time)
  c = as.vector(i$con)
  cl = as.vector(i$conlog)
  tr = as.vector(i$treatm)
  plot(t, c, pch = 19, col=ifelse(i$treat == "REF", "blue", "red"), grid(10), xlab = "Time after dosing (h)", ylab = "Plasma Drug concentration (pg/mL)",
       main = paste("Subject ", i$sub[1]))
  lines(t,c,type ="l")
  legend("topright", legend = c("Reference", "Test"), col = c("blue","red"), pch = c(19,17), bty = "o")
  plot(t, cl, pch = 19, col=ifelse(i$treat == "REF", "blue", "red"), grid(10), xlab = "Time after dosing (h)", ylab = "Log(Plasma Drug concentration (pg/mL))",
       main = paste("Log(Subject ", i$sub[1],")"))
  lines(t,cl,type ="l")
  legend("topright", legend = c("Reference", "Test"), col = c("blue","red"), pch = c(19,17), bty = "o")
}

## b. SPAGHETTI PLOT
r = subset(table, treatm == "REF")
t = subset(table, treatm == "TEST")
## b.1 Reference
plot(r$time, r$con, type = "l", col = "blue", grid(10),
       xlab = "Time after dosing (h)", ylab = "Plasma Reference concentration (pg/mL)",
       main = paste("Spaghetti plot - Reference"))
## b.2 Logreference
plot(r$time, log(r$con), type = "l", col = "blue", grid(10),
     xlab = "Time after dosing (h)", ylab = "Log10(Plasma Reference concentration (pg/mL))",
     main = paste("Spaghetti plot - Log10(Reference)"))

## b.3 Test
plot(t$time, t$con, type = "l", col = "red", grid(10),
     xlab = "Time after dosing (h)", ylab = "Plasma Test concentration (pg/mL)",
     main = paste("Spaghetti plot - Test"))
## b.4 Logtest
plot(t$time, log(t$con), type = "l", col = "red", grid(10),
     xlab = "Time after dosing (h)", ylab = "Log10(Plasma Test concentration (pg/mL))",
     main = paste("Spaghetti plot - Log10(Test)"))

## c. MEAN CONCENTRATION PLOTS

me.sp = data.frame()
for (i in 1:nrow(sp$`1`)){
  for (j in sp){
   me.sp = rbind(me.sp,j[i])
  }
}

me.sp2 = split(me.sp, me.sp$theo)

meansp = data.frame()
for (i in me.sp2){
  meansp = c(meansp, split(i, i$treatm))
}

add.ref = vector()
add.tes = vector()
add.ref.log = vector()
add.tes.log = vector()
for (j in meansp){
  if (j$treatm[1] == "REF"){
    add.ref = c(add.ref, sum(j$con, na.rm = T)/nrow(j))
    add.ref.log = c(add.ref.log, sum(j$conlog, na.rm = T)/nrow(j))
  } else {
    add.tes = c(add.tes, sum(j$con, na.rm = T)/nrow(j))
    add.tes.log = c(add.tes.log, sum(j$conlog, na.rm = T)/nrow(j))
  }
}

theo2 = unique(theo)
## LINEAR
plot(theo2, add.ref, pch = 19, col = "blue", grid(10),
     xlab = "Time post-dose (h)", ylab = "Mean plasma Drug Concentration(pg/mL)",
     main = paste("Drug Mean plasma concentration"))
lines(theo2, add.ref, col = "blue")
par(new = T)
plot(theo2, add.tes, pch = 17 , col = "red", grid(10), xlab = "", ylab = "", axes = F)
lines(theo2, add.tes, col = "red")
legend("topright", legend = c("Reference", "Test"), col = c("blue","red"), pch = c(19,17), bty = "o")
## SEMILOGARITHMIC
plot(theo2, add.ref.log, pch = 19, col = "blue", grid(10),
     xlab = "Time post-dose (h)" , ylab = "Log(Mean plasma Drug Concentration(pg/mL))" ,
     main = paste("Log Drug Mean Plasma Concentration"))
lines(theo2, add.ref.log, col = "blue")
par(new = T)
plot(theo2, add.tes.log, pch = 17 , col = "red", grid(10), xlab = "" , ylab = "", axes = F )
lines(theo2, add.tes.log, col = "red")
legend("topright", legend = c("Reference", "Test"), col = c("blue","red"), pch = c(19,17), bty = "o")

dev.off()
## We open a PDF which will store the results of the loop. We transform t and c into vectors (no se porqué, pero sino no funciona)
## In the plot the dots will be blue if the treatment is "REF" and red if it is the "TEST"
## We close the device outside the loop so it will create only one pdf and append all the graphs at the end

### PREPARE SUMMARY STATISTICS DATA FRAME
stat = rbind("Cmax_Ref", "Cmax_Test", "AUC_Ref", "AUC_Test")

### 3. CALCULATE CMAX FOR EACH SUBJECT AND TREATMENT ###
mc = data.frame()
sp = split(table, list(table$treatm, table$subj))

tmax = vector()
for (i in sp){
    cmax = max(i$con, na.rm = T)
    mc = unique(rbind(mc, data.frame(i$subj, i$seq, i$treatm, cmax)))
}

colnames(mc) = c("Subject","Sequence","Treatment","Cmax(pg/mL)")

### STATISTICS

## Metric + Nr. of subjects
nr.sub = nrow(mc)/2

### 3.a ARITHMETIC MEAN CMAX for REF and TEST

mean.cmax.ref = mean(mc$`Cmax(pg/mL)`[mc$Treatment == "REF"])
mean.cmax.tes = mean(mc$`Cmax(pg/mL)`[mc$Treatment == "TEST"])

### 3.b SD CMAX for REF and TEST

sd.cmax.ref = sd(mc$`Cmax(pg/mL)`[mc$Treatment == "REF"])
sd.cmax.tes = sd(mc$`Cmax(pg/mL)`[mc$Treatment == "TEST"])

## 3.c MEDIAN CMAX for REF and TEST

med.cmax.ref = median(mc$`Cmax(pg/mL)`[mc$Treatment == "REF"])
med.cmax.tes = median(mc$`Cmax(pg/mL)`[mc$Treatment == "TEST"])

## 3.d GEOMETRIC MEAN CMAX for REF and TEST
geo.cmax.ref = exp(mean(log(mc$`Cmax(pg/mL)`[mc$Treatment == "REF"])))
geo.cmax.tes = exp(mean(log(mc$`Cmax(pg/mL)`[mc$Treatment == "TEST"])))

## 3.e % CV CMAX for REF and TEST (SD/MEAN)*100
cv.cmax.ref = (sd.cmax.ref/mean.cmax.ref)*100
cv.cmax.tes = (sd.cmax.tes/mean.cmax.tes)*100

## 3.f MIN CMAX for REF and TEST
min.cmax.ref = min(mc$`Cmax(pg/mL)`[mc$Treatment == "REF"])
min.cmax.tes = min(mc$`Cmax(pg/mL)`[mc$Treatment == "TEST"])

## 3.g MAX CMAX for REF and TEST
max.cmax.ref = max(mc$`Cmax(pg/mL)`[mc$Treatment == "REF"])
max.cmax.tes = max(mc$`Cmax(pg/mL)`[mc$Treatment == "TEST"])

## ORGANISE RESULTS TOGETHER
cmax.ref.res = cbind("Cmax_Ref",nr.sub, mean.cmax.ref, sd.cmax.ref, med.cmax.ref, geo.cmax.ref, cv.cmax.ref, min.cmax.ref, max.cmax.ref )
cmax.tes.res = cbind("Cmax_Test",nr.sub, mean.cmax.tes, sd.cmax.tes, med.cmax.tes, geo.cmax.tes, cv.cmax.tes, min.cmax.tes, max.cmax.tes)

### 4. AUC FOR EACH PATIENT AND TREATMENT ###

## a) AUC with the linear trapezoidal method is calculated through this formula:
## AUC = sum (0.5 * (t2-t1) * (C1+C2))

## a.1) Add consecutive concentration values (for all the DFs)
## Output is one numeric vector with the results of all the DFs 
co = vector()
resco = numeric()
for (i in sp){
    co = as.vector(i$con)
    for (j in 1:length(co)-1){
      resco = c(resco, sum(co[j] + co[j+1]))
    }
}

## Add the sumed concentrations to the DF and repeat the splitting command to update sp with the new column
table$conadd = resco
sp = split(table, list(table$treatm, table$subj))

## a.2) Substract consecutive time values from the general table. In order not to substract
## time values from different subjects, an NA is introduced when the value is negative

resti = 0
tim = table$time
for (i in 1:(length(tim)-1)){
  resti = c(resti, tim[i+1]-tim[i])
  if (resti[i] < 0){
    resti[i] = NA
  }
}

## Add the substracted times to the DF and repeat the splitting command to update sp with the new column
table$timesub = resti
sp = split(table, list(table$treatm, table$subj))

## a.3) Calculate the area of each trapezoid (0.5 * (t2-t1) * (C1+C2)) 

timesub = table$timesub
conadd = table$conadd
trp = vector()

for (i in 1:length(timesub)){
    trp = c(trp, 0.5*timesub[i]*conadd[i])
}

## Add a column with the area of each trapezoid to the table DF
table$lintrapezoid = trp[1:length(trp)]
sp = split(table, list(table$treatm, table$subj))

## a.4) AUC for each subject and treatment
auc = vector()
for (i in sp){
  auc = c(auc, sum(i$lintrapezoid,na.rm = T))
}

## Attach auc to mc which already contains the values for cmax (mc will have the Primary PK Parameters)

mc = cbind(mc, auc)
colnames(mc) = c("Subject","Sequence","Treatment","Cmax(pg/mL)","AUClast")

write.table(mc, "PrimaryPKParameters.txt", row.names = FALSE, sep="\t", 
            dec=",")

### STATISTICS

## a.5) Mean AUC for REF and TEST

mean.auc.ref = mean(mc$AUClast[mc$Treatment == "REF"])
mean.auc.tes = mean(mc$AUClast[mc$Treatment == "TEST"])

## a.6) SD AUC for REF and TEST

sd.auc.ref = sd(mc$AUClast[mc$Treatment == "REF"])
sd.auc.tes = sd(mc$AUClast[mc$Treatment == "TEST"])

## a.7) MEDIAN AUC for REF and TEST
med.auc.ref = median(mc$AUClast[mc$Treatment == "REF"])
med.auc.tes = median(mc$AUClast[mc$Treatment == "TEST"])

## a.8) GEOMETRIC MEAN AUC for REF and TEST
geo.auc.ref = exp(mean(log(mc$AUClast[mc$Treatment == "REF"])))
geo.auc.tes = exp(mean(log(mc$AUClast[mc$Treatment == "TEST"])))

## a.9) % CV AUC for REF and TEST (SD/MEAN)*100
cv.auc.ref = (sd.auc.ref/mean.auc.ref)*100
cv.auc.tes = (sd.auc.tes/mean.auc.tes)*100

## a.10) MIN AUC for REF and TEST
min.auc.ref = min(mc$AUClast[mc$Treatment == "REF"])
min.auc.tes = min(mc$AUClast[mc$Treatment == "TEST"])

## a.11) MAX AUC for REF and TEST
max.auc.ref = max(mc$AUClast[mc$Treatment == "REF"])
max.auc.tes = max(mc$AUClast[mc$Treatment == "TEST"])

## ORGANISE RESULTS TOGETHER
auc.ref.res = cbind("AUC_Ref", nr.sub, mean.auc.ref, sd.auc.ref, med.auc.ref, geo.auc.ref, cv.auc.ref, min.auc.ref, max.auc.ref )
auc.tes.res = cbind("AUC_Test", nr.sub, mean.auc.tes, sd.auc.tes, med.auc.tes, geo.auc.tes, cv.auc.tes, min.auc.tes, max.auc.tes )

## b) AUC with the logarithmic trapezoidal method is calculated through this formula:
## AUC = sum ((C1+C2)/(ln(C1)-ln(C2)) * (t2-t1))
## log in R calculates by default natural logarithms (log neperiano)

## b.1) Transform to logarithm all the concentration values

clog = log(table$con)
table$conlog = clog

## b.2) Substract consecutive concentration logarithmic values (C1-C2)

sublog = vector()
conlog = table$conlog
for (i in 1:length(conlog)){
  sublog = c(sublog, conlog[i]-conlog[i+1])
}

table$sublog = sublog

## b.3) Substract consecutive concentration values (C1-C2)
consub = vector()
for (i in 1:length(con)){
  consub = c(consub, con[i]-con[i+1])
}

table$consub = consub

## b.3) Calculate the area of each trapezoid (C1-C2)/(ln(C1)-ln(C2)) * (t2-t1))

trpz = vector()
for (i in 1:length(sublog)){
  trpz = c(trpz, ((consub[i]/sublog[i])*timesub[i]))
}

table$logtrapezoid = trpz

## b.4) Sum the areas of all the trapezoids for each subject and treatment
sp = split(table, list(table$treatm, table$subj))
auclog = vector()
for (i in sp){
  auclog = c(auclog, sum(i$logtrapezoid,na.rm = T))
}

## c) AUC with linear-logarithmic trapezoidal method 
# Linear for the increasing []s and log for the decreasing []

linlog = vector()
auc.linlog = vector()
for (j in sp){
  for (i in 1:length(j$con)){
      if (is.na(con[i]== T)){
        next()
      }
      if (con[i] > con[i+1] || is.na(con[i+1] == T)){
        linlog = c(linlog, j$logtrapezoid[i])
      } else if (con[i] < con[i+1]){
        linlog = c(linlog, j$lintrapezoid[i])
      }
    
  }
  print(linlog)
  add = sum(linlog, na.rm = TRUE)
  auc.linlog = c(auc.linlog, add)
  linlog = vector()
}

##Alternative
cmax2 = mc$`Cmax(pg/mL)`
linlog = vector()
for (i in sp){
 for (j in i$con){
  if (cmax2 == j ){
    linlog = c(linlog, i$lintrapezoid[1:j])
  }
 }
}




### SUMMARY STATISTICS
stat = rbind(cmax.ref.res, cmax.tes.res, auc.ref.res, auc.tes.res)
colnames(stat) = c("Metric_Treatm","n","ArithmMean", "SD", "Median", "GeomMean", "CV(%)","Min","Max")

## GENERATE OUTPUT FILE

write.table(stat, "summary_stats.txt", row.names = F, sep="\t", dec=",")

### 5. ANOVA FOR CMAX AND AUC BETWEEN REF AND TEST ###
## a. ANOVA for CMAX

library(nlme) ##package for non-mixed effects
cmax.log = log(mc$Cmax)
mc = cbind(mc, cmax.log) 

## Linear mixed effects model including Sequence nested into Subject as a random effect, 
## the rest have fixed effects
cmax.mod = lme(cmax.log ~ Subject + Sequence + Treatment, 
              random = ~1 | Subject/Sequence, data = mc, na.action = na.exclude)

summary(cmax.mod)
### 90% CONFIDENCE INTERVALS
cmax.int = intervals(cmax.mod, level = 0.9, which = "fixed")

cmax.df = as.data.frame(c(cmax.int))
cmax.lower = exp(cmax.df[4,1])*100
cmax.upper = exp(cmax.df[4,3])*100

cmax.res = as.data.frame(cbind(cmax.lower, cmax.upper))

cmax.an = anova(cmax.mod)

## a1. IS BIOEQUIVALENT?
if (cmax.lower > 80 && cmax.upper < 125){
  cmax.res = cbind(cmax.res, "YES")
} else{
  cmax.res = cbind(cmax.res, "NO")
}

## b. ANOVA for AUC
library(nlme) ##package for non-mixed effects
auc.log = log(mc$AUClast)
mc$AUC.log = auc.log

auc.mod = lme(AUC.log ~ Subject + Sequence + Treatment, 
               random = ~1 | Subject/Sequence, data = mc, na.action = na.exclude)

summary(auc.mod)
auc.an = anova(auc.mod)
auc.int = intervals(auc.mod, level = 0.9, which = "fixed")

auc.df = as.data.frame(c(auc.int))
auc.lower = exp(auc.df[4,1])*100
auc.upper = exp(auc.df[4,3])*100

auc.res = as.data.frame(cbind(auc.lower, auc.upper))

### b1. IS BIOEQUIVALENT?
if (auc.lower > 80 && auc.upper < 125){
  auc.res = cbind(auc.res, "YES")
} else{
  auc.res = cbind(auc.res, "NO")
}
  
### 6.RATIO CALCULATION FOR LN-TRANSFORMED DATA
## RATIO(%REF) = 100exp(TestLSM - RefLSM)

## a. CMAX RATIO
library(lsmeans)
lsm.cmax = lsmeans(cmax.mod, ~Treatment)
s = summary(lsm.cmax)

lsm.cmax = exp(mean(log(cmax.mod, ~Treatment)))
## Obtain LSMean values for REF and TEST and convert them to numeric
lsm.cmax.s = as.data.frame(s[c("lsmean")]) 
var = as.numeric(lsm.cmax.s[,1])

## Calculate ratio
r.cmax = 100*exp(var[2]-var[1]) 

cmax.res = cbind(r.cmax, cmax.res )
cmax.res = cbind("Ln(Cmax)", cmax.res)
colnames(cmax.res) = c("Dependent","%Ratio", "Lower_CI_90", "Upper_CI_90", "Bioequivalence")

## b. AUC RATIO
library(lsmeans)
lsm.auc = lsmeans(auc.mod, ~Treatment)
su = summary(lsm.auc)

## b. AUC RATIO (geometric means)

## Obtain LSMean values for REF and TEST and convert them to numeric
lsm.auc.s = as.data.frame(su[c("lsmean")]) 
var1 = as.numeric(lsm.auc.s[,1])

## Calculate ratio
r.auc = 100*exp(var1[2]-var1[1])

auc.res = cbind(r.auc, auc.res)
auc.res = cbind("Ln(AUC)", auc.res)
colnames(auc.res) = c("Dependent","%Ratio", "Lower_CI_90", "Upper_CI_90", "Bioequivalence")

### 7. GENERATE RESULTS FILE
### a. FILE
results = rbind(cmax.res, auc.res)

write.table(results, "bioequivalence_res.txt", 
            row.names = FALSE, sep="\t", dec=",")
