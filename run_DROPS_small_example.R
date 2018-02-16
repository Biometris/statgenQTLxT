

setwd('M:/willem/research/STATISTICAL_GENETICS/multi_env_gwas_paper/code/git')

source('estimateEffects.R')

source('make.V.inv.array.R')

source('LL.quad.form.diag.R')

load(file='DROPS_small_example.RData')

r1 <- estimateEffects(W = W, mm = mm, Y = Y, Vg = Vg, Ve = Ve, K = K)


