library("edgeR")
source("/ufrc/barbazuk/lboat/Old_World_New_BED/trago_output/ase_counts/parents/rnaseq_plot_funcs.R")

setwd('/ufrc/barbazuk/lboat/Old_World_New_BED/trago_output/ase_counts/parents/Tpr_Tdu_voom')
targets <- readTargets()
print(targets)

d <- readDGE(targets, sep=",")
print(dim(d))
#print(head(d))
colnames(d) <- c("Tdu1","Tdu2","Tdu3","Tpr1","Tpr2","Tpr3","Tms1","Tms2","Tms3")

all_samples = as.matrix(d$counts)
colnames(all_samples) <- d$samples$description
write.csv(all_samples, "all_samples.csv")

keep <- rowSums(cpm(d) > cpm(10,mean(d$samples$lib.size))[1]) >= 9
d <- d[keep,]

pdf("additive_boxplot_log-CPM.pdf")
boxplot(cpm(d$counts + 1, log=TRUE))
dev.off()

#d <- calcNormFactors(d, method = "TMM")

trt = factor(c(1,1,1,2,2,2,3,3,3), labels=c("Tdu","Tpr","Tms"))
design=model.matrix(~0 + trt)

pdf("additive_voom_plot.pdf")
voom=voom(d, design, plot=TRUE)
dev.off()

write.table(voom$E, "additive_voom_expression_values.txt", sep="\t", quote=F, row.names = TRUE)
fit = lmFit(voom, design)
overall_model <- eBayes(fit)

pdf("additive_residual_std_dev.pdf")
plotSA(overall_model)
dev.off()

topTable(overall_model, coef=ncol(design))
top=topTable(overall_model, sort="none", n=Inf, coef=ncol(design))
write.table(top, "additive_DE_overall_model.txt", sep="\t", quote=F, row.names = TRUE)

print(summary(top))

# Results are relative to Tdu
contrast.matrix = makeContrasts(contrasts=c("trtTpr-trtTdu"), levels=design)
fit2 = contrasts.fit(fit, contrast.matrix)
fit2 = eBayes(fit2)

topTable(fit2, coef=ncol(contrast.matrix))
top=topTable(fit2, sort="none", n=Inf, coef=ncol(contrast.matrix))
write.table(top, "DE_Tdu_Tpr.txt", sep="\t", quote=F, row.names = TRUE)

contrast.matrix=makeContrasts(contrasts=c("(trtTdu+trtTpr)/2-trtTms"), levels=design)
fit2=contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

#qqt(fit2$t , df=fit2$df.prior+fit2$df.residual,pch=16,cex=0.2,main="(Tdu+Tpo)/2 - Tm Q-Q plot")
#abline(0,1)

topTable(fit2, coef=ncol(contrast.matrix))
top=topTable(fit2, sort="none", n=Inf, coef=ncol(contrast.matrix))

# plot_MA_and_Volcano(logCounts = top$AveExpr,logFoldChange = top$logFC,FDR = top$adj.P.Val,xlab="logCounts",ylab = "logFC", title = "MA plot")

write.table(top, "DE_additive_Tdu_Tpr-Tms.txt", sep="\t", quote=F, row.names = TRUE)

pdf("additive_MDS_plot.pdf")
plotMDS(voom, main="MDS plot")
dev.off()
