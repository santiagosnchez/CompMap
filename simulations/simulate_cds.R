# scale : approximate median size
get_gene_sizes <- function(N, scale=1000, offset=300){
    # random integers from gamma distribution
    gene_len = round(rgamma(n=N, scale=1000, shape=1.35), 0)
    gene_len = gene_len - (gene_len %% 3) # make them multples of 3
    while (sum(gene_len < 300) != 0){
        repacement = round(rgamma(n=sum(gene_len < 300), scale=1000, shape=1.35), 0)
        repacement = repacement = repacement - (repacement %% 3)
        gene_len[ which(gene_len < 300) ] = repacement
    }
    return(gene_len)
}

# all possible codons
make_codons <- function(){
    codons = list()
    nucleotides = c("T","C","A","G")
    for (first in nucleotides)
        for (second in nucleotides)
            for (third in nucleotides)
                codons$codons = append(codons$codons,paste(first,second,third,sep=""))
    codons$start = "ATG"
    codons$stop = codons$codons[grep("TAA|TAG|TGA", codons$codons)]
    codons$codons = codons$codons[-grep("TAA|TAG|TGA", codons$codons)]
    codons$fourfold = codons$codons[grep("[ATGC]C[ATGC]|[CG]T[ATGC]|GT[ATGC]|[CG]G[ATGC]", codons$codons)]
    codons$twofold = codons$codons[grep("TT[ATGC]|TA[TC]|CA[ATGC]|AA[ATGC]|GA[ATGC]|TG[TC]|AG[ATGC]", codons$codons)]
    codons$threefold = codons$codons[grep("AT[ATC]", codons$codons)]
    codons$synonymous = c(codons$fourfold, codons$twofold, codons$threefold)
    return(codons)
}

# create a synonymous mutational model
synonymous_mutations <- function(codons){
    mutations = list()
    # start with two-fold
    for (cod in codons$twofold){
        if (substr(cod, 3, 3) == "T"){
            mutations[[cod]] = sub("T$","C", cod)
        } else if (substr(cod, 3, 3) == "C"){
            mutations[[cod]] = sub("C$","T", cod)
        } else if (substr(cod, 3, 3) == "G"){
            mutations[[cod]] = sub("G$","A", cod)
        } else if (substr(cod, 3, 3) == "A"){
            mutations[[cod]] = sub("A$","G", cod)
        }
    }
    # then with four-fold
    for (cod in codons$fourfold){
        if (substr(cod, 3, 3) == "T"){
            for (j in c("A","C","G")){
                mutations[[cod]] = c(mutations[[cod]], sub("T$", j, cod))
            }
        } else if (substr(cod, 3, 3) == "C"){
            for (j in c("A","T","G")){
                mutations[[cod]] = c(mutations[[cod]], sub("C$", j, cod))
            }
        } else if (substr(cod, 3, 3) == "G"){
            for (j in c("A","C","T")){
                mutations[[cod]] = c(mutations[[cod]], sub("G$", j, cod))
            }
        } else if (substr(cod, 3, 3) == "A"){
            for (j in c("T","C","G")){
                mutations[[cod]] = c(mutations[[cod]], sub("A$", j, cod))
            }
        }
    }
    # three-fold last
    for (cod in codons$threefold){
        if (substr(cod, 3, 3) == "T"){
            for (j in c("A","C")){
                mutations[[cod]] = c(mutations[[cod]], sub("T$", j, cod))
            }
        } else if (substr(cod, 3, 3) == "C"){
            for (j in c("A","T")){
                mutations[[cod]] = c(mutations[[cod]], sub("C$", j, cod))
            }
        } else if (substr(cod, 3, 3) == "A"){
            for (j in c("T","C")){
                mutations[[cod]] = c(mutations[[cod]], sub("A$", j, cod))
            }
        }
    }
    return(mutations)
}

# sample codons to make CDS
make_cds <- function(codons, gene_sizes){
    sequences = lapply(gene_sizes, function(x, cod=codons){
        # codons are randomly sampled, adding start and stop codons, no premature stops
        cds = paste(c(cod$start, sample(cod$codons, x-2, replace=T), sample(cod$stop,1)), collapse="")
        return(cds)
    })
    return(sequences)
}

seq_to_codons <- function(cds){
    substring(cds, seq(1, nchar(cds), 3), seq(3, nchar(cds), 3))
}

# substitute synonymous codons given a synonymous substitution rate
synonymous_divergence <- function(cds, codons, mutations, rate){
    cds_1 = seq_to_codons(cds)
    cds_2 = cds_1
    total_syn = length(cds_1[ cds_1 %in% codons$synonymous ])
    nsites = round(total_syn * rate, 0)
    codon_positions = which(cds_1 %in% codons$synonymous)
    if (nsites < total_syn){
        codon_positions = sample(codon_positions, nsites, replace=F)
    }
    for (i in codon_positions){
        cds_2[i] = sample(mutations[[cds_1[i]]], 1)
    }
    return(paste(cds_2, collapse=""))
}

# simulate data
# N = number of genes
# scale = approx. mean of the protein length in aa
# offset minimum protein size
# CDS is a list
simulate_synonymous_divergence <- function(N, basename, scale=1000, offset=300, rate=0.01){
    codons = make_codons()
    mutations = synonymous_mutations(codons)
    aa_gene_sizes = get_gene_sizes(N, scale, offset)
    CDS_1 = make_cds(codons, aa_gene_sizes)
    CDS_2 = lapply(CDS_1, synonymous_divergence, codons=codons, mutations=mutations, rate=rate)
    ii = 0
    cat("",file=paste(basename,N,rate,"d1.fasta",sep="_"))
    for (i in CDS_1){
        ii = ii + 1
        cat(c(paste(paste0(">d1.g",ii),"t1",sep="."),i), sep="\n", file=paste(basename,N,rate,"d1.fasta",sep="_"), append=T)
    }
    ii = 0
    cat("",file=paste(basename,N,rate,"d2.fasta",sep="_"))
    for (i in CDS_2){
        ii = ii + 1
        cat(c(paste(paste0(">d2.g",ii),"t1",sep="."),i), sep="\n", file=paste(basename,N,rate,"d2.fasta",sep="_"), append=T)
    }
}
