library("edgeR")
source("/ufrc/barbazuk/lboat/Old_World_New_BED/trago_output/ase_counts/parents/rnaseq_plot_funcs.R")

setwd('/ufrc/barbazuk/lboat/Old_World_New_BED/trago_output/ase_counts/parents/Tpr_Tdu_voom')
targets <- readTargets()
print(targets)

d <- readDGE(targets, sep=",")
print(dim(d))
#print(head(d))
colnames(d) <- c("Tdu_1","Tdu_2","Tdu_3","Tpr_1","Tpr_2","Tpr_3")

keep <- rowSums(cpm(d) > cpm(10,mean(d$samples$lib.size))[1]) >= 6
d <- d[keep,]

pdf("boxplot_log-CPM.pdf")
boxplot(cpm(d$counts + 1, log=TRUE))
dev.off()

#d <- calcNormFactors(d, method = "TMM")

trt = factor(c(1,1,1,2,2,2), labels=c("Tdu","Tpr"))
design=model.matrix(~0 + trt)
voom=voom(d, design, plot=TRUE)

write.table(voom$E, "voom_expression_values.txt", sep="\t", quote=F, row.names = TRUE)
fit = lmFit(voom, design)
overall_model <- eBayes(fit)

pdf("residual_std_dev.pdf")
plotSA(overall_model)
dev.off()

topTable(overall_model, coef=ncol(design))
top=topTable(overall_model, sort="none", n=Inf, coef=ncol(design))
write.table(top, "DE_overall_model.txt", sep="\t", quote=F, row.names = TRUE)

print(summary(top))

contrast.matrix = makeContrasts(contrasts=c("trtTdu-trtTpr"), levels=design)
fit2 = contrasts.fit(fit, contrast.matrix)
fit2 = eBayes(fit2)
#qqt(fit2$t , df=fit2$df.prior+fit2$df.residual,pch=16,cex=0.2,main="Tdu-Tpo Q-Q plot")
#abline(0,1)
topTable(fit2, coef=ncol(contrast.matrix))
top=topTable(fit2, sort="none", n=Inf, coef=ncol(contrast.matrix))
write.table(top, "DE_Tdu_Tpr.txt", sep="\t", quote=F, row.names = TRUE)

#contrast.matrix = makeContrasts(contrasts=c("trtTpo-trtTdu"), levels=design)
#fit2 = contrasts.fit(fit, contrast.matrix)
#fit2 <- eBayes(fit2)

#dt = decideTests(fit2)
#print(dt)
#write.fit(fit2, dt, "DE_Tpo-Tdu.txt", adjust="BH")
#topTable(fit2, coef=ncol(design))
#top = topTable(fit2, sort="none", n=Inf, coef=ncol(design))
#write.table(top, "DE_Tpo-Tdu.txt", sep="\t", quote=F, row.names = TRUE)

