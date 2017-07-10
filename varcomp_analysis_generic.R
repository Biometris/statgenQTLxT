######################################
# START OF THE ANALYSIS


# re-run model without any SNP, for comparison of variance components
# Important: we use 'fixed.formula' defined above

#setwd(wd2)

if (use.pc) {

  pc.part <- paste0('PC',1:n.pc)

  pc.part <- paste(pc.part,collapse=' + ')

  fixed.formula <- paste0(trait.name,' ~ Env + ', pc.part , ' + Env:(',pc.part,')')

  fixed.formula.red <- paste0(trait.name,' ~ Env + ', pc.part)

} else {

  fixed.formula <- paste0(trait.name,' ~ Env + S + Env:S')

  fixed.formula.red <- paste0(trait.name,' ~ Env + S')

}

snp.part <- paste(g$snps$snp,collapse=' + ')

####################################
# Model (1) from the draft

fixed.formula1 <- fixed.formula

reml.obj1 <- asreml(fixed  = as.formula(fixed.formula1),
                    random = ~ G + G:L + G:Y + G:L:Y,
                    rcov   = ~ at(Env):units,
                    data   = d,
                    maxit  = 500,
                    workspace   = workspace.size)

wald(reml.obj1)

####################################
# Model (2) from the draft

fixed.formula2 <- fixed.formula

reml.obj2 <- asreml(fixed  = as.formula(fixed.formula2),
                    random = ~ G + G:EC + G:L + G:Y + G:L:Y,
                    rcov   = ~ at(Env):units,
                    data   = d,
                    maxit  = 500,
                    workspace   = workspace.size)

####################################
# Model (3) from the draft

fixed.formula3 <- paste0(fixed.formula,' + Env:(',snp.part,')')

reml.obj3 <- asreml(fixed  = as.formula(fixed.formula3),
                    random = ~ G + G:EC + G:L + G:Y + G:L:Y,
                    rcov   = ~ at(Env):units,
                    data   = d,
                    maxit  = 500,
                    workspace   = workspace.size)

#reml.obj3.pred <- predict(object=reml.obj3, classify='G',
#                                 workspace   = workspace.size, data   = d)

MU <- mean(reml.obj3$fitted)

# A model with G.all indices; didn't work very well
#reml.obj2b <- asreml(fixed  = as.formula(fixed.formula3),
#                    random = ~ G + G:(indice.ET0+ indice.Light+ indice.MaxTemp+ indice.AvTemp+ indice.DayTemp+ indice.NightTemp+ indice.MaxVPD+ indice.nbHsup30deg+ indice.SoilPsi) + G:L + G:Y + G:L:Y,
#                    rcov   = ~ at(Env):units,
#                    data   = d,
#                    maxit  = 50,
#                    workspace   = workspace.size)

indices <- c('indice.ET0','indice.Light','indice.MaxTemp','indice.AvTemp','indice.DayTemp',
             'indice.NightTemp','indice.MaxVPD','indice.nbHsup30deg','indice.SoilPsi','treatment')

n.indices <- length(indices)

treatment <- rep('irrigated',nrow(d))
treatment[grep(d$Env,pattern='drought')] <- 'drought'
d$treatment <- factor(treatment)
# table(d$treatment,d$hyd)

##############################
# Indices (one each time), without SNPs and without G.EC

m2b       <- list()
m2b.table <- matrix(NA,5,n.indices)

for (j in 1:n.indices) {

  ft <- call('as.formula',paste0(' ~G + G:', indices[j],' + G:L + G:Y + G:L:Y'))

  reml.obj2b <- asreml(fixed  = as.formula(fixed.formula2),
                       random = eval(ft),
                       rcov   = ~ at(Env):units,
                       data   = d,
                       maxit  = 50,
                       workspace   = workspace.size)

  m2b[[j]]        <- summary(reml.obj2b)$varcomp[1:5,2]
  names(m2b[[j]]) <- unlist(lapply(strsplit(rownames(summary(reml.obj2b)$varcomp[1:5,]),split='!'),function(x){x[1]}))

  m2b.table[,j]   <- as.numeric(summary(reml.obj2b)$varcomp[1:5,2])
}

colnames(m2b.table) <- indices
rownames(m2b.table) <- unlist(lapply(strsplit(rownames(summary(reml.obj2b)$varcomp[1:5,]),split='!'),function(x){x[1]}))

####################################
# Indices (one each time), without SNPs and with G.EC

m2c       <- list()
m2c.table <- matrix(NA,6,n.indices)

for (j in 1:n.indices) {

  ft <- call('as.formula',paste0(' ~G + G:EC + G:', indices[j],' + G:L + G:Y + G:L:Y'))

  reml.obj2c <- asreml(fixed  = as.formula(fixed.formula2),
                       random = eval(ft),
                       rcov   = ~ at(Env):units,
                       data   = d,
                       maxit  = 50,
                       workspace   = workspace.size)

  m2c[[j]]        <- summary(reml.obj2c)$varcomp[1:6,2]
  names(m2c[[j]]) <- unlist(lapply(strsplit(rownames(summary(reml.obj2c)$varcomp[1:6,]),split='!'),function(x){x[1]}))

  m2c.table[,j]   <- as.numeric(summary(reml.obj2c)$varcomp[1:6,2])
}

colnames(m2c.table) <- indices
rownames(m2c.table) <- unlist(lapply(strsplit(rownames(summary(reml.obj2c)$varcomp[1:6,]),split='!'),function(x){x[1]}))



##############################
# Indices (one each time), with SNPs and without G.EC

m3b       <- list()
m3b.table <- matrix(NA,5,n.indices)

for (j in 1:n.indices) {

  ft <- call('as.formula',paste0(' ~G + G:', indices[j],' + G:L + G:Y + G:L:Y'))

  reml.obj3b <- asreml(fixed  = as.formula(fixed.formula3),
                       random = eval(ft),
                       rcov   = ~ at(Env):units,
                       data   = d,
                       maxit  = 50,
                       workspace   = workspace.size)

  m3b[[j]]        <- summary(reml.obj3b)$varcomp[1:5,2]
  names(m3b[[j]]) <- unlist(lapply(strsplit(rownames(summary(reml.obj3b)$varcomp[1:5,]),split='!'),function(x){x[1]}))

  m3b.table[,j]   <- as.numeric(summary(reml.obj3b)$varcomp[1:5,2])
}

colnames(m3b.table) <- indices
rownames(m3b.table) <- unlist(lapply(strsplit(rownames(summary(reml.obj3b)$varcomp[1:5,]),split='!'),function(x){x[1]}))

####################################
# Indices (one each time), with SNPs and with G.EC

m3c       <- list()
m3c.table <- matrix(NA,6,n.indices)

for (j in 1:n.indices) {

  ft <- call('as.formula',paste0(' ~G + G:EC + G:', indices[j],' + G:L + G:Y + G:L:Y'))

  reml.obj3c <- asreml(fixed  = as.formula(fixed.formula3),
                       random = eval(ft),
                       rcov   = ~ at(Env):units,
                       data   = d,
                       maxit  = 50,
                       workspace   = workspace.size)

  m3c[[j]]        <- summary(reml.obj3c)$varcomp[1:6,2]
  names(m3c[[j]]) <- unlist(lapply(strsplit(rownames(summary(reml.obj3c)$varcomp[1:6,]),split='!'),function(x){x[1]}))

  m3c.table[,j]   <- as.numeric(summary(reml.obj3c)$varcomp[1:6,2])
}

colnames(m3c.table) <- indices
rownames(m3c.table) <- unlist(lapply(strsplit(rownames(summary(reml.obj3c)$varcomp[1:6,]),split='!'),function(x){x[1]}))

rownames(m3b.table)[2] <- rownames(m2b.table)[2] <- rownames(m3c.table)[3] <- rownames(m2c.table)[3] <- 'G.index'



####################################
# Model (4) from the draft

fixed.formula4 <- paste0(fixed.formula.red,' + Env:(',snp.part,')')

reml.obj4 <- asreml(fixed  = as.formula(fixed.formula4),
                    random = ~ G + G:EC + G:L + G:Y + G:L:Y,
                    rcov   = ~ at(Env):units,
                    data   = d,
                    maxit  = 500,
                    workspace   = workspace.size)

####################################


# without SNPs, with E:PC, without G.EC
m1 <- summary(reml.obj1)$varcomp[1:4,2]
names(m1) <- unlist(lapply(strsplit(rownames(summary(reml.obj1)$varcomp[1:4,]),split='!'),function(x){x[1]}))

# without SNPs, with E:PC, with G.EC
m2 <- summary(reml.obj2)$varcomp[1:5,2]
names(m2) <- unlist(lapply(strsplit(rownames(summary(reml.obj2)$varcomp[1:5,]),split='!'),function(x){x[1]}))

# with SNPs, with E:PC, with G.EC
m3 <- summary(reml.obj3)$varcomp[1:5,2]
names(m3) <- unlist(lapply(strsplit(rownames(summary(reml.obj3)$varcomp[1:5,]),split='!'),function(x){x[1]}))

# with SNPs, with E:PC, with G.EC
#m2b <- summary(reml.obj2b)$varcomp[1:13,2]
#names(m2b) <- unlist(lapply(strsplit(rownames(summary(reml.obj2b)$varcomp[1:13,]),split='!'),function(x){x[1]}))

# with SNPs, without E:PC, with G.EC
m4 <- summary(reml.obj4)$varcomp[1:5,2]
names(m4) <- unlist(lapply(strsplit(rownames(summary(reml.obj4)$varcomp[1:5,]),split='!'),function(x){x[1]}))

#############
# print output

# capture.output(source('varcomp_analysis_generic.R'),file=output.file)

#zz <- file(output.file, open = "wt")

#sink(zz)

sink(file=output.file) # ,type='message')

cat('=====================================================','\n\n')
cat('TRAIT: ',trait.name,'\n\n','Estimated variance components','\n\n')
cat('=====================================================','\n\n')



cat('Model 1','\n\n')

print(m1)

cat('\n\n')



cat('Model 2 (without SNPs, with G.EC)','\n\n')

print(m2)

cat('\n\n')



cat('Model 2b (Indices, without SNPs and without G.EC)','\n\n')

print(m2b.table)

cat('\n\n')




cat('Model 2c (Indices, without SNPs and with G.EC)','\n\n')

print(m2c.table)

cat('\n\n')



cat('Model 3 (with SNPs and with G.EC)','\n\n')

print(m3)

cat('\n\n')


cat('Model 3b (Indices, with SNPs and without G.EC)','\n\n')

print(m3b.table)

cat('\n\n')



cat('Model 3c (Indices, with SNPs and with G.EC)','\n\n')

print(m3c.table)

cat('\n\n')






cat('Model 4 (with SNPs and with G.EC, but WITHOUT Env:PC)','\n\n')

print(m4)

cat('\n\n')


cat('=====================================================','\n\n')
cat('TRAIT: ',trait.name,'\n\n','Square root of variance components','\n\n')
cat('=====================================================','\n\n')



cat('Model 1','\n\n')

print(sqrt(m1))

cat('\n\n')



cat('Model 2 (without SNPs, with G.EC)','\n\n')

print(sqrt(m2))

cat('\n\n')



cat('Model 2b (Indices, without SNPs and without G.EC)','\n\n')

print(sqrt(m2b.table))

cat('\n\n')



cat('Model 2c (Indices, without SNPs and with G.EC)','\n\n')

print(sqrt(m2c.table))

cat('\n\n')





cat('Model 3 (with SNPs and with G.EC)','\n\n')

print(sqrt(m3))

cat('\n\n')



cat('Model 3b (Indices, with SNPs and without G.EC)','\n\n')

print(sqrt(m3b.table))

cat('\n\n')


cat('Model 3c (Indices, with SNPs and with G.EC)','\n\n')

print(sqrt(m3c.table))

cat('\n\n')




cat('Model 4 (with SNPs and with G.EC, but WITHOUT Env:PC)','\n\n')

print(sqrt(m4))

cat('\n\n')


cat('=====================================================','\n\n')
cat('TRAIT: ',trait.name,'\n\n','Square root of variance components, divided by the mean, times 100 percent','\n\n')
cat('=====================================================','\n\n')

cat('\n\n', 'MEAN: ', MU,'\n\n\n')

cat('Model 1','\n\n')

print(100 * sqrt(m1) / MU )

cat('\n\n')

cat('Model 2 (without SNPs, with G.EC)','\n\n')

print(100 * sqrt(m2) / MU )

cat('\n\n')



cat('Model 2b (Indices, without SNPs and without G.EC)','\n\n')

print(100 * sqrt(m2b.table) / MU )

cat('\n\n')



cat('Model 2c (Indices, without SNPs and with G.EC)','\n\n')

print(100 * sqrt(m2c.table) / MU )

cat('\n\n')



cat('Model 3 (with SNPs and with G.EC)','\n\n')

print(100 * sqrt(m3) / MU )

cat('\n\n')

cat('Model 3b (Indices, with SNPs and without G.EC)','\n\n')

print(100 * sqrt(m3b.table) / MU )

cat('\n\n')




cat('Model 3c (Indices, with SNPs and with G.EC)','\n\n')

print(100 * sqrt(m3c.table) / MU )

cat('\n\n')


#cat('Model 2b','\n\n')
#print(100 * sqrt(m2b) / MU )
#cat('\n\n')

cat('Model 4 (with SNPs and with G.EC, but WITHOUT Env:PC)','\n\n')

print(100 * sqrt(m4) / MU )

cat('\n\n')



sink(NULL)


###############################
