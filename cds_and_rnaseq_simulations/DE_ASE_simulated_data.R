# DE and stats
library(DESeq2)
library(limma)
library(edgeR)
library(MASS)
library(car)
# data wrangling
library(tidyr)
library(dplyr)
library(tibble)
# plotting
library(gridExtra)
library(RColorBrewer)
library(ggsci)
library(ggplot2)
library(cowplot)
library(lemon)
library(ggridges)
# multithread processing
library(parallel)
library(latex2exp)

theme_set(theme_cowplot())

# get file list
files = dir("counts/final", full.names=T)
sims = c("highrate","lowrate","mediumrate")
groups = c("species","alleles","compmap")

# read all data
all_data = list()
for (sim in sims){
    all_data[[sim]] = list()
    sim_files = files[ grep(sim, files) ]
    for (group in groups){
        all_data[[sim]][[group]] = read.csv(sim_files[ grep(group, sim_files) ], sep="\t", row.names=1)
    }
}
#str(all_data)

# tidy data.frames
gene_names = rownames(all_data[[sims[1]]][[groups[1]]])
for (sim in sims){
    for (group in groups){
        # same row order
        tmp = all_data[[sim]][[group]]
        colnames(tmp) = gsub("04","01", colnames(tmp))
        colnames(tmp) = gsub("05","02", colnames(tmp))
        colnames(tmp) = gsub("06","03", colnames(tmp))
        all_data[[sim]][[group]] = tmp[gene_names,]
        # reorder columns in compmap output
        if (group == "compmap"){
            tmp = all_data[[sim]][[group]]
            tmp = tmp[,c(c(1,3,5),c(2,4,6))]
            colnames(tmp) = paste("compmap", colnames(tmp), sep="_")
            all_data[[sim]][[group]] = tmp
        }
    }
}
rm(tmp)
#str(all_data)

# make coldata table
sample_names = as.vector(sapply(all_data[[sims[1]]], colnames))
coldata = data.frame(samples=sample_names, group=substr(sample_names, 1, 7), dataset=substr(sample_names, 19, 20))

# plot all counts
biplots = list()
for (sim in sims){
    tmp1 = all_data[[sim]][["alleles"]]
    tmp2 = all_data[[sim]][["compmap"]]
    compare_counts = cbind(tmp1,tmp2, gene=rownames(tmp1))
    df.compare_counts = compare_counts %>%
        pivot_longer(cols=-gene, names_to="samples", values_to="counts") %>%
        extract(samples, c("group","sample","dataset"), "(alleles|compmap)_(sample_0[1-6])_(d1|d2)") %>%
        pivot_wider(names_from=group, values_from=counts) %>%
        mutate(ratio = log2(compmap) / log2(alleles))
    biplots[[sim]] = ggplot(df.compare_counts, aes(log2(alleles), log2(compmap), color=sample)) +
        geom_point(show.legend=F) +
        geom_smooth(method="lm", se=F, color="black", linetype=2) +
        facet_rep_grid(sample~dataset) +
        labs(x="True log2 counts (featureCounts)",y="Competitive mapping log2 counts (CompMap)") +
        background_grid()
}
biplots[["highrate"]] = biplots[["highrate"]] + theme(axis.title.x=element_blank()) + ggtitle("high rate (0.1 changes/site)")
biplots[["mediumrate"]] = biplots[["mediumrate"]] + theme(axis.title.y=element_blank()) + ggtitle("medium rate (0.01 changes/site)")
biplots[["lowrate"]] = biplots[["lowrate"]] + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + ggtitle("low rate (0.001 changes/site)")
plot_grid(biplots[["highrate"]], biplots[["mediumrate"]], biplots[["lowrate"]], nrow=1, align="h")
ggsave("figures/compare_log2_counts.png")

# consolidate
all_data_merge = data.frame()
for (sim in sims){
    tmp1 = all_data[[sim]][["alleles"]]
    tmp2 = all_data[[sim]][["compmap"]]
    colnames(tmp1) = colnames(tmp2) = sub("alleles_|compmap_","", colnames(tmp1))
    tmp1$genes = tmp2$genes = rownames(tmp1)
    tmp1$sim = tmp2$sim = sim
    tmp1$source = "alleles"
    tmp2$source = "compmap"
    all_data_merge = rbind(all_data_merge, tmp1, tmp2)
}
all_data_merge = all_data_merge %>% pivot_longer(cols=c(-genes, -sim, -source), values_to="counts", names_to="sample")
all_data_merge = all_data_merge %>% extract(col=sample, into=c("rep","dataset"), regex=c("sample_(0[1-3])_(d[12])"))
one_sample = all_data_merge %>% filter(rep == "01" & dataset == "d1")
one_sample = one_sample %>% pivot_wider(names_from=source, values_from=counts)
# rate_names = c(expression(paste(italic("d")[s],"=0.1")), expression(paste(italic("d")[s],"=0.01")), expression(paste(italic("d")[s],"=0.001")))
# names(rate_names) = c("highrate","mediumrate","lowrate")
#one_sample$sim = factor(c("high div. (ps = 0.1)", "low rate (0.001)", "medium rate (0.01)")[ factor(one_sample$sim) ], c("low rate (0.001)", "medium rate (0.01)", "high rate (0.1)"))
one_sample$sim = factor(one_sample$sim, c("lowrate","mediumrate","highrate"))
one_sample$sim = factor(c("0.001","0.01","0.1")[one_sample$sim], c("0.001","0.01","0.1"))
# one_sample$sim = rate_names[one_sample$sim]

appender <- function(string)
    TeX(paste("$\\textit{\\d_{\\s}}=$", string))

counts_comp = ggplot(one_sample, aes(log2(compmap), log2(alleles))) +
    geom_point(color="grey40", alpha=0.1) +
    geom_abline(slope=1, color="black", linetype=1, size=1) +
    scale_x_continuous(limits=c(6, 16), expand=c(0,0)) +
    scale_y_continuous(limits=c(6, 16), expand=c(0,0)) +
    labs(x=expression(paste("log"[2]," read counts [CompMap]")),
         y=expression(atop(paste("log"[2]," read counts"),"[featureCounts]"))) +
    facet_rep_wrap(~sim, labeller = as_labeller(appender, default = label_parsed)) +
    background_grid()

# differential expression
species = alleles = compmap = vector("list", 3)
names(species) = names(alleles) = names(compmap) = sims
for (sim in sims){
    # species
    dat = all_data[[sim]][["species"]]
    colnames(dat) = sub("species_","", colnames(dat))
    dds = DESeqDataSetFromMatrix(countData = dat, colData = coldata[ coldata$group == "species", ], design = ~ dataset)
    dds = DESeq(dds)
    species[[sim]] = results(dds)
    # alleles
    dat = all_data[[sim]][["alleles"]]
    colnames(dat) = sub("alleles_","", colnames(dat))
    dds = DESeqDataSetFromMatrix(countData = dat, colData = coldata[ coldata$group == "alleles", ], design = ~ dataset)
    dds = DESeq(dds)
    alleles[[sim]] = results(dds)
    # compmap
    dat = all_data[[sim]][["compmap"]]
    colnames(dat) = sub("compmap_","", colnames(dat))
    dds = DESeqDataSetFromMatrix(countData = dat, colData = coldata[ coldata$group == "compmap", ], design = ~ dataset)
    dds = DESeq(dds)
    compmap[[sim]] = results(dds)
    #
    # species[[sim]] = all_data[[sim]][["species"]]
    # colnames(species[[sim]]) = sub("species_","", colnames(species[[sim]]))
    # alleles[[sim]] = all_data[[sim]][["alleles"]]
    # colnames(alleles[[sim]]) = sub("alleles_","", colnames(alleles[[sim]]))
    # compmap[[sim]] = all_data[[sim]][["compmap"]]
    # colnames(compmap[[sim]]) = sub("compmap_","", colnames(compmap[[sim]]))
}

# difference in logfold change

diff_lfc = data.frame(diff=c(alleles[[sims[[1]]]]$log2FoldChange-compmap[[sims[[1]]]]$log2FoldChange,
                            alleles[[sims[[2]]]]$log2FoldChange-compmap[[sims[[2]]]]$log2FoldChange,
                            alleles[[sims[[3]]]]$log2FoldChange-compmap[[sims[[3]]]]$log2FoldChange),
                     sim=rep(sims, each=1000))
diff_lfc$sim = factor(c("high rate (0.1)", "low rate (0.001)", "medium rate (0.01)")[ factor(diff_lfc$sim) ], c("low rate (0.001)", "medium rate (0.01)", "high rate (0.1)"))

hist_diff = ggplot(diff_lfc, aes(diff)) +
    geom_histogram(fill="black") +
    facet_rep_wrap(~sim) +
    labs(x=expression(paste("log"[2],"-FC [featureCounts] - log"[2],"-FC [CompMap]")), y="count (genes)") +
    facet_rep_wrap(~sim) +
    background_grid()

error_rate = data.frame(diff=c((alleles[[sims[[1]]]]$padj < 0.05 & compmap[[sims[[1]]]]$padj < 0.05) | (alleles[[sims[[1]]]]$padj > 0.05 & compmap[[sims[[1]]]]$padj > 0.05),
                            (alleles[[sims[[2]]]]$padj < 0.05 & compmap[[sims[[2]]]]$padj < 0.05) | (alleles[[sims[[2]]]]$padj > 0.05 & compmap[[sims[[2]]]]$padj > 0.05),
                            (alleles[[sims[[3]]]]$padj < 0.05 & compmap[[sims[[3]]]]$padj < 0.05) | (alleles[[sims[[3]]]]$padj > 0.05 & compmap[[sims[[3]]]]$padj > 0.05)),
                     sim=rep(sims, each=1000))
false_positives = data.frame(diff=c(alleles[[sims[[1]]]]$padj > 0.05 & compmap[[sims[[1]]]]$padj < 0.05,
                            alleles[[sims[[2]]]]$padj > 0.05 & compmap[[sims[[2]]]]$padj < 0.05,
                            alleles[[sims[[3]]]]$padj > 0.05 & compmap[[sims[[3]]]]$padj < 0.05),
                     sim=rep(sims, each=1000))
false_negatives = data.frame(diff=c(alleles[[sims[[1]]]]$padj < 0.05 & compmap[[sims[[1]]]]$padj > 0.05,
                            alleles[[sims[[2]]]]$padj < 0.05 & compmap[[sims[[2]]]]$padj > 0.05,
                            alleles[[sims[[3]]]]$padj < 0.05 & compmap[[sims[[3]]]]$padj > 0.05),
                     sim=rep(sims, each=1000))
calc_error_rate = error_rate %>% group_by(sim, diff) %>% count() %>% na.omit() %>% group_by(sim) %>% summarize(total=(n/sum(n))*100) %>% filter(total < 20)
calc_error_rate = bind_cols(calc_error_rate, false_positives %>% group_by(sim, diff) %>%
                                                count() %>% na.omit() %>% group_by(sim) %>%
                                                summarize(false_positive=(n/sum(n))*100) %>% filter(false_positive < 20) %>%
                                                ungroup() %>% select(false_positive))
calc_error_rate = bind_cols(calc_error_rate, false_negatives %>% group_by(sim, diff) %>%
                                                count() %>% na.omit() %>% group_by(sim) %>%
                                                summarize(false_negative=(n/sum(n))*100) %>% filter(false_negative < 20) %>%
                                                ungroup() %>% select(false_negative))
calc_error_rate = calc_error_rate %>% pivot_longer(cols=-sim, names_to="type_error", values_to="error") %>% mutate(type_error=sub("_"," ",type_error))
calc_error_rate = calc_error_rate %>% mutate(sim=factor(sim, c("lowrate","mediumrate","highrate")))

overall_error = ggplot(calc_error_rate, aes(x=sim, shape=type_error, group=type_error)) +
    geom_point(aes(y=error), size=3, position=position_dodge(0.5), show.legend=F) +
    #geom_text(aes(label=error, y=error+0.1), position=position_dodge(0.5)) +
    #geom_linerange(aes(y=error, ymin=0, ymax=error), position=position_dodge(0.5)) +
    scale_y_log10(breaks=c(0.01,0.1,1,10), limits=c(0.01,50), labels=c("0.01","0.1","1","10"), expand=c(0,0)) +
    scale_shape_manual(values=15:17, name="") +
    scale_x_discrete(label=c("0.001", "0.01", "0.1")) +
    labs(x=TeX("$\\textit{\\d_{\\s}}$"), y=TeX("error rate ($\\log_{\\10}$)")) + #        expression(atop(paste("error rate"),paste("(% log"[10], " scale)")))) +
    annotation_logticks(sides="l") +
    background_grid(major="y") +
    theme(legend.position=c(.6,.9), legend.text=element_text(size=10))

# cis and trans classification

# test for differences in log2 fold change differences
# based on a Wald-type posthoc test (linearHypothesis) and
# a linear model
getTransEffects = function(x, reps=3){
	x = as.data.frame(x)
	x$T = rep(c("P","H"), each=reps*2)
	x$S = rep(rep(c("sp1","sp2"), each=reps), 2)
	fit = lm(x ~ T/S, data=x)
	fit2 = car::linearHypothesis(fit, c("TH:Ssp2 = TP:Ssp2"))
	p.value = fit2$`Pr(>F)`[2]
	list(x, fit, fit2, p.value)
}

# trans: vector of adjusted pvalues
# allele: ASE results from DESeq2
# species: differential expression analyses between species from DESeq2
classify = function(trans, allele, species){
    x = rep(NA, 1000)
    x[ species$padj > 0.05 & allele$padj > 0.05 ] = "conserved"
    x[ is.na(x) & trans > 0.05 & allele$padj < 0.05 & species$padj < 0.05 ] = "cis only"
    x[ is.na(x) & trans < 0.05 & allele$padj > 0.05 & species$padj < 0.05 ] = "trans only"
    x[ is.na(x) & trans < 0.05 & allele$padj < 0.05 & species$padj > 0.05 ] = "compensatory"
    x[ is.na(x) & trans < 0.05 & ((allele$log2FoldChange > 0 & species$log2FoldChange < 0) | (allele$log2FoldChange < 0 & species$log2FoldChange > 0)) ] = "cis x trans"
    x[ is.na(x) & trans < 0.05 & ((allele$log2FoldChange > 0 & species$log2FoldChange > 0) | (allele$log2FoldChange < 0 & species$log2FoldChange < 0)) ] = "cis + trans"
    x[ is.na(x) ] = "ambiguous"
    return(x)
}

# get trans values

trans.fc = trans.cm = vector("list", 3)
names(trans.fc) = names(trans.cm) = sims
for (sim in sims){
    dat = cpm(cbind(all_data[[sim]][["species"]], all_data[[sim]][["alleles"]]), log=T)
    print("running fc...")
    res = mclapply(split(dat, 1:dim(dat)[1]), getTransEffects, mc.cores = 4)
    pval = sapply(res, function(x) x[[4]] )
    pval = p.adjust(pval, method="BH")
    trans.fc[[sim]] = pval
    dat = cpm(cbind(all_data[[sim]][["species"]], all_data[[sim]][["compmap"]]), log=T)
    print("running cp...")
    res = mclapply(split(dat, 1:dim(dat)[1]), getTransEffects, mc.cores = 4)
    pval = sapply(res, function(x) x[[4]] )
    pval = p.adjust(pval, method="BH")
    trans.cm[[sim]] = pval
}

# classify

classification.fc = classification.cm = vector("list", 3)
names(classification.fc) = names(classification.cm) = sims
for (sim in sims){
    classification.fc[[sim]] = classify(trans.fc[[sim]], alleles[[sim]], species[[sim]])
    classification.cm[[sim]] = classify(trans.cm[[sim]], compmap[[sim]], species[[sim]])
}
classes = sort(unique(classification.fc[[1]])[-3])

# consolidate error

error_by_class = data.frame()
for (sim in sims){
    df = data.frame(sim=sim)
    x = data.frame(fc=classification.fc[[sim]], cm=classification.cm[[sim]])
    df = cbind(df, as.data.frame(x %>% group_by(fc) %>% summarise(false_negative=1-(sum(fc == cm)/length(fc))) %>% ungroup() %>% filter( fc != "ambiguous")))
    df = cbind(df, as.data.frame(x %>% group_by(cm) %>% summarise(false_positive=1-(sum(cm == fc)/length(cm))) %>% ungroup() %>% filter( cm != "ambiguous")))
    error_by_class = rbind(error_by_class, df[,-4])
}
error_by_class = error_by_class %>% pivot_longer(cols=c(-sim, -fc), names_to="error_type", values_to="error") %>% rename("class" = fc)
error_by_class = error_by_class %>% mutate(sim=factor(sim, c("lowrate","mediumrate","highrate")))

# plot
mean_error_per_class = ggplot(error_by_class %>% filter(error_type=="false_negative"), aes(x=sim)) +
    geom_point(aes(y=error*100), position=position_jitter(width=0.1), alpha=0.5) +
    geom_errorbar(data = error_by_class %>% filter(error_type=="false_negative") %>% group_by(sim) %>% summarize(mean_error=mean(error)),
        aes(ymax=mean_error * 100, ymin=mean_error * 100), width=0.4, size=1.5) +
    #geom_linerange(aes(ymax=error * 100, ymin=0), position=position_dodge(0.5)) +
    scale_y_log10(breaks=c(0.01,0.1,1,10), limits=c(0.01,50), labels=c("0","0.1","1","10"), expand=c(0,0)) +
    scale_x_discrete(label=c("0.001", "0.01", "0.1")) +
    labs(x=TeX("$\\textit{\\d_{\\s}}$"), y=TeX("error rate ($\\log_{\\10}$)")) + #
    annotation_logticks(sides="l") +
    #annotate(geom="text", label="zero", x=2.5, y=0.05) +
    #annotate(geom="segment", x=2.45, y=0.03, xend=2.1, yend=0.01, arrow = arrow(length = unit(0.2, "cm"))) +
    #annotate(geom="segment", x=2.55, y=0.03, xend=2.9, yend=0.01, arrow = arrow(length = unit(0.2, "cm"))) +
    annotation_logticks(sides="l") +
    background_grid(major="y")
    #theme(axis.text.x=element_text(angle=90, size=9, hjust=1, vjust=0.5))

# gline = linesGrob(y = c(0.013, 0.015), x = c(-.015, .015),  gp = gpar(col = "black", lwd = 2.5))
# gline2 = linesGrob(y = c(-0.25, 0.5),x = c(0, 0),  gp = gpar(col = "red", lwd = 5))


top = plot_grid(counts_comp +
                    theme(strip.background=element_blank()),
                hist_diff +
                    theme(strip.background=element_blank(), strip.text=element_blank()), ncol=1, align="v", labels=c("a","b"))

bottom = plot_grid(overall_error, mean_error_per_class + theme(axis.title.y=element_blank()), labels=c("c", ""))
plot_grid(top, bottom, rel_heights=c(0.6,0.3), ncol=1)

ggsave(file="figure.pdf", device="pdf", useDingbats=F)
ggsave(file="figure_alt.png", dpi=300)
