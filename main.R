library("vcfR")
library("ggplot2")

cosmicMu_coding <- read.vcfR("CosmicCodingMuts.vcf")
cosmicMu_coding <- cosmicMu_coding@fix
cosmicMu_coding <- as.data.frame(cosmicMu_coding, stringsAsFactors = F)
tableRes <- data.frame(table(cosmicMu_coding$CHROM))

# factor(tableRes$Var1, levels = tableRes$Var1[sort(tableRes$Freq, decreasing = T, index.return = T)[[2]]])
tableRes$Var1 <- factor(tableRes$Var1, levels = tableRes$Var1[order(tableRes$Freq, decreasing = T)])
ggplot(data = tableRes, aes(x = tableRes$Var1, y = tableRes$Freq, 
                            fill = tableRes$Var1)) + geom_bar(stat = "identity")

tableRes$Freq[24] / sum(tableRes$Freq[1:22]) * 100

cosmicMu_ChrX_coding <- cosmicMu_coding[cosmicMu_coding$CHROM == "X",]

rm(tableRes, cosmicMu_coding)
gc()

selectChrX <- function(inputPath) {
        cosmicMu <- read.vcfR(inputPath)
        cosmicMu <- cosmicMu@fix
        cosmicMu <- as.data.frame(cosmicMu, stringsAsFactors = F)
        tableRes <- data.frame(table(cosmicMu$CHROM))
        
        tableRes$Var1 <- factor(tableRes$Var1, levels = tableRes$Var1[order(tableRes$Freq, decreasing = T)])
        ggplot(data = tableRes, aes(x = tableRes$Var1, y = tableRes$Freq, 
                                    fill = tableRes$Var1)) + geom_bar(stat = "identity")
        
        message("Chromosome X Percentage: ")
        cat(tableRes$Freq[24] / sum(tableRes$Freq[1:22]) * 100)
        
        cosmicMu_ChrX <- cosmicMu[cosmicMu$CHROM == "X",]
        cosmicMu_ChrX
}

cosmicMu_ChrX_noncoding <- selectChrX("CosmicNonCodingVariants.vcf")

tmpData <- read.table("1000GenomesProj_SNV_ChrX.tsv", sep = "\t", stringsAsFactors = F)

processData <- function(inputRow) {
        inputALT <- inputRow[[4]]
        inputALT_vec <- strsplit(inputALT, ",")[[1]]
        copyRow <- inputRow[rep(1, length(inputALT_vec)),]
        copyRow$ALT <- inputALT_vec
        copyRow
}

total <- 15055
pb <- txtProgressBar(min = 0, max = total, style = 3)
for(row.idx in 1:total){
        inputRow <- multiData[row.idx,]
        outputRow <- processData(inputRow)
        if (row.idx == 1) output_multi <- outputRow else output_multi <- rbind(output_multi, outputRow)
        setTxtProgressBar(pb, row.idx)
}

ncharData <- sapply(tmpData$ALT, nchar, USE.NAMES = F)
uniData <- tmpData[ncharData == 1,]

finalData <- rbind(uniData, output_multi)
finalData_2 <- finalData[sort(finalData$POS, index.return = T)$ix,]

cosmicMu_ChrX <- cosmicMu_ChrX[sort(cosmicMu_ChrX$POS, index.return = T)$ix,]

selectREF <- sapply(cosmicMu_ChrX$REF, nchar, USE.NAMES = F)
table(selectREF)
cosmicMu_ChrX_2 <- cosmicMu_ChrX[selectREF == 1,]

selectALT <- sapply(cosmicMu_ChrX_2$ALT, nchar, USE.NAMES = F)
table(selectALT)
cosmicMu_ChrX_2 <- cosmicMu_ChrX_2[selectALT == 1,]

cosmicMu_ChrX_2$POS <- as.character(cosmicMu_ChrX_2$POS)
GenomeProj_processed$POS <- as.character(GenomeProj_processed$POS)

ref_cosmic <- apply(cosmicMu_ChrX_2[,2:4], 1, paste0, collapse = ".")
alt_genome <- apply(GenomeProj_processed[,2:4], 1, paste0, collapse = ".")

overlap_idx <- alt_genome %in% ref_cosmic
table(overlap_idx)

GenomeProj_processed_2 <- GenomeProj_processed[!overlap_idx,]

chk_genome <- apply(GenomeProj_processed_2[,2:4], 1, paste0, collapse = ".")

check_idx <- chk_genome %in% ref_cosmic
table(check_idx)

save(GenomeProj_processed_2, file = "GenomeProj_processed.RData")


cosmicMu_ChrX_2$POS <- as.numeric(cosmicMu_ChrX_2$POS)
GenomeProj_processed_2$POS <- as.numeric(GenomeProj_processed_2$POS)

altPos <- cosmicMu_ChrX_2$POS
startIDX <- altPos[1]
rangeIDX <- seq(startIDX - 1000, startIDX + 1000, 1)
pd <- txtProgressBar(min = 2, max = length(altPos), style = 3)
for (i in 2:length(altPos)) {
        if(i %% 100 == 0) message(i)
        newRangeIDX <- seq(altPos[i] - 1000, altPos[i] + 1000, 1)
        rangeIDX <- unique(c(rangeIDX, newRangeIDX))
        # setTxtProgressBar(pd, i)
}

outbound <- 1
rangeIDX <- seq(altPos[1] - 1000, altPos[1] + 1000, 1)
for (part in 1:1034){
        start <- outbound + 1
        outbound <- part * 1000
        message(start, ":", outbound)
        if(outbound > 1033816) outbound <- 1033816
        for (i in start:outbound) {
                message("  ", i)
                newRangeIDX <- seq(altPos[i] - 1000, altPos[i] + 1000, 1)
                rangeIDX <- unique(c(rangeIDX, newRangeIDX))
        }
        
}

rangeIDX <- seq(altPos[1] - 1000, altPos[1] + 1000, 1)
res <- sapply(2:length(altPos), function(x) {
        if(x %% 100 == 0) message(x)
        updateIDX <- get("rangeIDX", envir = environment())
        newRangeIDX <- seq(altPos[i] - 1000, altPos[i] + 1000, 1)
        rangeIDX <<- unique(c(updateIDX, newRangeIDX))
        # setTxtProgressBar(pd, i)
})

which(cosmicMu_ChrX_2$POS > cosmicMu_ChrX_2$POS[1] + 1000)

infoPattern = c("AC", "AF", "AFR_AF", "AMR_AF", "AN", "DP", 
                "EAS_AF", "EUR_AF", "NS", "SAS_AF", "VT")

infoPattern = c("GENE", "STRAND", "CDS", "AA", "CNT")

INFO_sep <- sapply(1:total, function(idx, infoPattern = infoPattern) {
        # print(idx)
        x <- df_fix_INFO[[idx]]
        tmp_res <- strsplit(x, ";")[[1]]
        
        tmp_res_2_2 <- regexpr("=", tmp_res)
        tmp_res_2_2 <- ifelse(tmp_res_2_2 == -1, F,T)
        tmp_res <- tmp_res[tmp_res_2_2]
        
        tmp_res_2 <-sapply(tmp_res, function(y) {
                tmp_res_3 <- strsplit(y, "=")[[1]]
        }, USE.NAMES = F)
        setTxtProgressBar(pb, idx)
        idx_check <- tmp_res_2[1,] %in% infoPattern
        res <- tmp_res_2[2,idx_check]
        names(res) <- tmp_res_2[1,idx_check]
        res
}, USE.NAMES = F)

GenomeProjSNV_single_ChrX_2$POS <- as.numeric(GenomeProjSNV_single_ChrX_2$POS)

tab.pos <- table(GenomeProjSNV_single_ChrX_2$POS)
tab.pos[which(tab.pos != 1)]

which(GenomeProjSNV_single_ChrX_2$POS == 61671)
GenomeProjSNV_single_ChrX_2[134:135,]

hg19_phast46 <- read.table("chrX.phastCons46way.wigFix", header = T, sep = "\n")
hg19_phast46$POS <- seq(60546, length.out = ncol(hg19_phast46))
hg19_phast46 <- hg19_phast46[,]

######

coding_cosmic <- cosmicData_hg38[cosmicData_hg38$REGION == "Coding",]
noncoding_cosmic <- cosmicData_hg38[cosmicData_hg38$REGION != "Coding",]
rm(cosmicData_hg38)

cd_cnt <- data.frame(table(coding_cosmic$CNT), region = "Coding", stringsAsFactors = F)
nc_cnt <- data.frame(table(noncoding_cosmic$CNT), region = "NonCoding", stringsAsFactors = F)

CNT_data <- rbind(cbind(coding_cosmic$CNT, coding_cosmic$REGION), cbind(noncoding_cosmic$CNT, noncoding_cosmic$REGION))

library(ggplot2)

data_cnt <- as.data.frame(data_cnt)
sapply(data_cnt, class)
data_cnt$Freq <- as.numeric(data_cnt$Freq) 

save(data_cnt, file = "data_cnt.RData")
write.csv(data_cnt, file = "data_cnt.csv")
data_cnt<-read.csv("ggplot_cnt.csv")
data_cnt$RecurrenceNumber <- factor(data_cnt$RecurrenceNumber, levels = c(seq(1:19), ">=20"))

ggplot(data = data_cnt, aes(x = RecurrenceNumber, y = Percentage, group = Region, 
                                 colour = Region)) +
        # geom_line(size = 1, aes(linetype = Region)) + 
        geom_line(size = 1) +
        scale_y_sqrt(breaks = c(90, 70, 50, 30, 15, 5, 2, 0.4, 0.05)) +
        theme(legend.position = "bottom",
              axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12),
              axis.title.x = element_text(size = 12),
              axis.title.y = element_text(size = 12),
              axis.text.y  = element_text(size = 12),
              legend.title = element_text(size = 12),
              legend.text  = element_text(size = 12)) # 5 * 12 in
save.image("ggplot_cnt_data.RData")

#####

GenomeData_hg38 <- read.table("GenomeProjData_hg38_cleaned.csv", header = T, sep = ",", stringsAsFactors = F)
phast30 <- read.table("hg38_phast30.csv", header = T, sep = ",", stringsAsFactors = F)
phylo30 <- read.table("hg38_phylo30.csv", header = T, sep = ",", stringsAsFactors = F)
phast100 <- read.table("hg38_phast100.csv", header = T, sep = ",", stringsAsFactors = F)
phylo100 <- read.table("hg38_phylo100.csv", header = T, sep = ",", stringsAsFactors = F)


addField <- function(GenomeData = GenomeData_hg38, conservScore = "phast30") {
        conserData <- get(conservScore)
        
        label <- conserData[[1]]
        names(label) <- conserData[[2]]
        levels_label <- as.list(label)
        
        input <- as.factor(GenomeData$POS)
        levels(input) <- levels_label
        
        res <- as.numeric(as.character(input))
        
        sample_idx <- sample(nrow(GenomeData), 1000)
        chk_res <- sapply(sample_idx, function(idx) {
                if (as.character(conserData[[2]][which(conserData[[1]] == GenomeData$POS[idx])]) == as.character(res[idx])) tmp <- T else tmp <- F
                tmp
        })
        if (!all(chk_res)) {
                stop("ERROR!")
        }
        
        GenomeData[[conservScore]] <- res
        GenomeData
}

GenomeData_hg38 <- addField(conservScore = "phast30")
GenomeData_hg38 <- addField(conservScore = "phast100") 
GenomeData_hg38 <- addField(conservScore = "phylo30")
GenomeData_hg38 <- addField(conservScore = "phylo100")

GenomeData_hg38$POS <- as.numeric(GenomeData_hg38$POS)
GenomeData_hg38$AF <- as.numeric(GenomeData_hg38$AF)

sapply(GenomeData_hg38, class)

write.csv(GenomeData_hg38, file = "GenomeProjData_hg38_v2_addScore.csv")

save(GenomeData_hg38, file = "GenomeProjData_hg38_v2_addScore.RData")

CosmicData_hg38 <- read.table("cosmicData_hg38_v1_cleaned.csv", header = T, sep = ",", stringsAsFactors = F)
phast30 <- read.table("COSMIC_hg38_phast30.csv", header = T, sep = ",", stringsAsFactors = F)
phylo30 <- read.table("COSMIC_hg38_phylo30.csv", header = T, sep = ",", stringsAsFactors = F)
phast100 <- read.table("COSMIC_hg38_phast100.csv", header = T, sep = ",", stringsAsFactors = F)
phylo100 <- read.table("COSMIC_hg38_phylo100.csv", header = T, sep = ",", stringsAsFactors = F)

CosmicData_hg38 <- addField(GenomeData = CosmicData_hg38, conservScore = "phast30")
CosmicData_hg38 <- addField(GenomeData = CosmicData_hg38, conservScore = "phast100") 
CosmicData_hg38 <- addField(GenomeData = CosmicData_hg38, conservScore = "phylo30")
CosmicData_hg38 <- addField(GenomeData = CosmicData_hg38, conservScore = "phylo100")

write.csv(CosmicData_hg38, file = "CosmicData_hg38_v2_addScore.csv")
save(CosmicData_hg38, file = "CosmicData_hg38_v2_addScore.RData")
save(GenomeData_hg38, file = "GenomeData_hg38_v2_addScore.RData")

POS_nc <- cosmicData_hg38$POS[cosmicData_hg38$REGION != "Coding"]

pb <- txtProgressBar(min = 0, max = 5301, style = 3)

POS_range <- seq(POS_nc[1] - 1000, POS_nc[1] + 1000)
sapply(2:5301, function(x) {
        setTxtProgressBar(pb, x)
        POS <- POS_nc[[x]]
        POS_range <<- (c(POS_range, seq(POS - 1000, POS + 1000)))
        if (x %% 500 == 0) POS_range <<- unique(POS_range)
        return(0)
})

GenomeData_hg38_nc <- GenomeData_hg38_window[GenomeData_hg38_window$REGION == "NonCoding",]
table(GenomeData_hg38_nc$POS %in% POS_range)

GenomeData_hg38_window$chk_POS <- "KEEP"
GenomeData_hg38_window$chk_POS[(GenomeData_hg38_window$REGION == "NonCoding") & (!GenomeData_hg38_window$POS %in% POS_range)] <- "OMIT"
table(GenomeData_hg38_window$chk_POS)

paste(GenomeData_hg38_window, , sep = ".") == row.names(genome_vep_aa_noncoding)

cosmicMu <- CosmicData_hg38[CosmicData_hg38$POS %in% names(which(table(CosmicData_hg38$POS)>1)),]
cosmicMu_2 <- CosmicData_hg38[CosmicData_hg38$POS %in% names(which(table(CosmicData_hg38$POS)==1)),]
dup_names <- names(which(table(CosmicData_hg38$POS)>1))

input_data <- cosmicMu
total <- length(dup_names)
# pb <- txtProgressBar(min = 0, max = total, style = 3)
refined_df <- sapply(1:length(dup_names), function(idx) {
        # setTxtProgressBar(pb, idx)
        x <- dup_names[[idx]]
        input_data_1 <- input_data[input_data$POS %in% x,]
        out_res <- c()
        
        for(ALT in unique(input_data_1$ALT)) {
                input_data_2 <- input_data_1[input_data_1$ALT %in% ALT,]
                res <- input_data_2[which.max(input_data_2$CNT),]
                out_res <- rbind(out_res, res)
        }
        out_res
})
refined_df_new <- apply(refined_df,1,unlist)
refined_df_new <- data.frame(refined_df_new, stringsAsFactors = F)
refined_df_new$POS <- as.numeric(refined_df_new$POS)
refined_df_new$CNT <- as.numeric(refined_df_new$CNT)


GenomeData_hg38_window <- GenomeData_hg38[GenomeData_hg38$POS %in% POS_range,]
save(GenomeData_hg38_window, file = "GenomeProjData_hg38_v3_filterWindow.RData")
write.csv(GenomeData_hg38_window, file = "GenomeProjData_hg38_v3_filterWindow.csv")


dataset <- rbind(cosmicData_hg38[,c(11, 7:10)], GenomeData_hg38_window[,c(10, 6:9)])


cosmic_vcf <- cbind(CHROM = cosmicData_hg38$CHROM, POS = cosmicData_hg38$POS, 
                    ID = paste(c(1:6561), cosmicData_hg38$label, cosmicData_hg38$REF, 
                               cosmicData_hg38$ALT, cosmicData_hg38$REGION, sep = "."), 
                    REF = cosmicData_hg38$REF, ALT = cosmicData_hg38$ALT, 
                    QUAL = ".", FILTER = ".", INFO = ".")
write.table(cosmic_vcf, file = "cosmic_hg38_v3.vcf", quote = F, sep = "\t", row.names = F, col.names = F)


genomeProj_vcf <- cbind(CHROM = GenomeData_hg38_window$CHROM, POS = GenomeData_hg38_window$POS, 
                        ID = paste(c(1:7700), GenomeData_hg38_window$label, GenomeData_hg38_window$REF, 
                                   GenomeData_hg38_window$ALT, sep = "."), 
                        REF = GenomeData_hg38_window$REF, ALT = GenomeData_hg38_window$ALT, 
                        QUAL = ".", FILTER = ".", INFO = ".")
write.table(genomeProj_vcf, file = "genomeProj_hg38_v3.vcf", quote = F, sep = "\t", row.names = F, col.names = F)

cosmic_vep <- read.table("vep_COSMIC", sep = "\t", header = T, stringsAsFactors = F)
cosmic_vep <- read.table("vep_CosmicCodingMuts_ChrX", sep = "\t", header = T, stringsAsFactors = F)

genome_vep <- read.table("vep_GenomeProj", sep = "\t", header = T, stringsAsFactors = F)


levels_consequence <- sapply(c(cosmic_vep$Consequence, genome_vep$Consequence), strsplit, split = ",")
levels_consequence <- unique(unlist(levels_consequence))


# input_data <- cosmic_vep

# AA_flag_name <- function(idx, input_data, AA.LETTERS = sort(c("G", "V", "L", "F", "P", "I", "A", "M", 
#                                                               "T", "S", "Y", "N", "Q", "W", "H", "K", 
#                                                               "R", "E", "D", "C"))) {
#         new_data <- input_data[input_data$Uploaded_variation == idx,]
#         res <- sapply(new_data$Amino_acids, AA_flag)
#         res <- apply(res, 1, max)
#         names(res) <- c(paste0("REF_", AA.LETTERS), paste0("ALT_", AA.LETTERS))
#         res
# }

vep_levels <- strsplit(("synonymous_variant, non_coding_transcript_exon_variant, missense_variant, stop_gained, frameshift_variant, incomplete_terminal_codon_variant, coding_sequence_variant, downstream_gene_variant, splice_region_variant, intron_variant, stop_retained_variant, upstream_gene_variant, splice_donor_variant, non_coding_transcript_variant, splice_acceptor_variant, 3_prime_UTR_variant, NMD_transcript_variant, 5_prime_UTR_variant, start_lost, inframe_deletion, stop_lost, inframe_insertion, protein_altering_variant, intergenic_variant, mature_miRNA_variant, start_retained_variant, transcript_ablation, transcript_amplification, TFBS_ablation, TFBS_amplification, TF_binding_site_variant, regulatory_region_ablation, regulatory_region_amplification, feature_elongation, regulatory_region_variant, feature_truncation"), 
                       ", ")[[1]]

AA_flag <- function(aminoAcidINFO, AA.LETTERS = sort(c("G", "V", "L", "F", "P", "I", "A", "M", 
                                                       "T", "S", "Y", "N", "Q", "W", "H", "K", 
                                                       "R", "E", "D", "C"))) {
        aa_vec <- strsplit(aminoAcidINFO, split = "/")[[1]]
        
        if (length(aa_vec) == 2) {
                if (aa_vec[[1]] == "*") REF_vec <- rep(1, 20) else REF_vec <- as.numeric(AA.LETTERS %in% aa_vec[[1]])
                if (aa_vec[[2]] == "*") ALT_vec <- rep(1, 20) else ALT_vec <- as.numeric(AA.LETTERS %in% aa_vec[[2]])
        } else if (aminoAcidINFO %in% AA.LETTERS) {
                if (aminoAcidINFO == "*") REF_vec <- rep(1, 20) else REF_vec <- as.numeric(AA.LETTERS %in% aminoAcidINFO[[1]])
                ALT_vec <- REF_vec
        } else {
                REF_vec <- rep(0, 20)
                ALT_vec <- REF_vec
        }
        out_vec <- c(REF_vec, ALT_vec)
        out_vec
}

all_flag_name <- function(idx_num, idx_name, input_data, vep_levels = vep_levels, 
                          AA.LETTERS = sort(c("G", "V", "L", "F", "P", "I", "A", "M", 
                                              "T", "S", "Y", "N", "Q", "W", "H", "K", 
                                              "R", "E", "D", "C"))) {
        idx <- idx_name[idx_num]
        if(idx_num %% 1000 == 0) message(idx_num, " / ", length(idx_name), " ",Sys.time())
        
        new_data <- input_data[input_data$Uploaded_variation == idx,]
        res <- sapply(new_data$Amino_acids, AA_flag, AA.LETTERS = AA.LETTERS)
        res <- apply(res, 1, max)
        # names(res) <- c(paste0("REF_", AA.LETTERS), paste0("ALT_", AA.LETTERS))
        aa_res <- res
        
        res <- lapply(new_data$Consequence, function(x) {
                out <- strsplit(x, ",")[[1]]
        })
        # browser()
        res <- unique(unlist(res))
        res <- table(factor(res, levels = vep_levels))
        cons_res <- c(num = nrow(new_data), res)
        
        output <- c(cons_res, aa_res)
        output
}
parallel::makeCluster(cl, 10)
parallel::clusterExport(cl, varlist = c("vep_levels", "AA_flag", "all_flag_name"))

cosmic_vep_all <- parallel::parSapply(cl, 1:length(unique(cosmic_vep$Uploaded_variation)), 
                         all_flag_name, idx_name = unique(cosmic_vep$Uploaded_variation), 
                         input_data = cosmic_vep, vep_levels = vep_levels,
                         AA.LETTERS = sort(c("G", "V", "L", "F", "P", "I", "A", "M", 
                                             "T", "S", "Y", "N", "Q", "W", "H", "K", 
                                             "R", "E", "D", "C")))
wrap_func <- function(idx_name, input_data = cosmic_vep) {
        
        vep_levels <- strsplit(("synonymous_variant, non_coding_transcript_exon_variant, missense_variant, stop_gained, frameshift_variant, incomplete_terminal_codon_variant, coding_sequence_variant, downstream_gene_variant, splice_region_variant, intron_variant, stop_retained_variant, upstream_gene_variant, splice_donor_variant, non_coding_transcript_variant, splice_acceptor_variant, 3_prime_UTR_variant, NMD_transcript_variant, 5_prime_UTR_variant, start_lost, inframe_deletion, stop_lost, inframe_insertion, protein_altering_variant, intergenic_variant, mature_miRNA_variant, start_retained_variant, transcript_ablation, transcript_amplification, TFBS_ablation, TFBS_amplification, TF_binding_site_variant, regulatory_region_ablation, regulatory_region_amplification, feature_elongation, regulatory_region_variant, feature_truncation"), 
                               ", ")[[1]]
        
        cosmic_vep_all <- sapply(1:length(idx_name), 
                                 all_flag_name, idx_name = idx_name, 
                                 input_data = input_data, vep_levels = vep_levels,
                                 AA.LETTERS = sort(c("G", "V", "L", "F", "P", "I", "A", "M", 
                                                     "T", "S", "Y", "N", "Q", "W", "H", "K", 
                                                     "R", "E", "D", "C")))
        cosmic_vep_all
}
cosmic_vep_all <- data.frame(t(cosmic_vep_all))
save(cosmic_vep_all, file = "cosmic_vep_all.RData")

# cons_fea <- function(idx, input_data, vep_levels = vep_levels) {
#         new_data <- input_data[input_data$Uploaded_variation == idx,]
#         res <- lapply(new_data$Consequence, function(x) {
#                 out <- strsplit(x, ",")[[1]]
#         })
#         # browser()
#         res <- unique(unlist(res))
#         res <- table(factor(res, levels = vep_levels))
#         output <- c(num = nrow(new_data), res)
#         return(output)
# }

cl <- parallel::makeCluster(28)
parallel::clusterExport(cl, varlist = "AA_flag")

cosmic_vep_aa <- parallel::parSapply(cl, unique(cosmic_vep$Uploaded_variation), AA_flag_name, input_data = cosmic_vep)
cosmic_vep_aa <- data.frame(t(cosmic_vep_aa))

genome_vep_aa <- sapply(unique(genome_vep$Uploaded_variation), AA_flag_name, input_data = genome_vep)
genome_vep_aa <- data.frame(t(genome_vep_aa))

save(cosmic_vep_aa, genome_vep_aa, file = "dataset_v1.2_vepAA.RData")

# test_vep <- read.table("vep_CosmicCodingMuts_ChrX", sep = "\t", header = T, stringsAsFactors = F)
# 
# vep_levels <- sapply(test_vep$Consequence, function(consequenceINFO) {
#         vep_cons <- strsplit(consequenceINFO, split = ",")[[1]]
#         vep_cons
# }, USE.NAMES = F)
# vep_levels <- c(unique(unlist(vep_levels)), "transcript_ablation", "transcript_amplification",
#                 "TFBS_ablation", "TFBS_amplification", "TF_binding_site_variant",
#                 "regulatory_region_ablation", "regulatory_region_amplification",
#                 "feature_elongation", "regulatory_region_variant", "feature_truncation")

vep_levels <- strsplit(("synonymous_variant, non_coding_transcript_exon_variant, missense_variant, stop_gained, frameshift_variant, incomplete_terminal_codon_variant, coding_sequence_variant, downstream_gene_variant, splice_region_variant, intron_variant, stop_retained_variant, upstream_gene_variant, splice_donor_variant, non_coding_transcript_variant, splice_acceptor_variant, 3_prime_UTR_variant, NMD_transcript_variant, 5_prime_UTR_variant, start_lost, inframe_deletion, stop_lost, inframe_insertion, protein_altering_variant, intergenic_variant, mature_miRNA_variant, start_retained_variant, transcript_ablation, transcript_amplification, TFBS_ablation, TFBS_amplification, TF_binding_site_variant, regulatory_region_ablation, regulatory_region_amplification, feature_elongation, regulatory_region_variant, feature_truncation"), 
                        ", ")[[1]]

tmp <- sapply(vep_levels, gsub, pattern = "_", replacement = " ", USE.NAMES = F)
write.csv(tmp, file = "vep_levels.csv", quote = F)

cons_fea <- function(idx, input_data, vep_levels = vep_levels) {
        new_data <- input_data[input_data$Uploaded_variation == idx,]
        res <- lapply(new_data$Consequence, function(x) {
                out <- strsplit(x, ",")[[1]]
        })
        # browser()
        res <- unique(unlist(res))
        res <- table(factor(res, levels = vep_levels))
        output <- c(num = nrow(new_data), res)
        return(output)
}

cosmic_vep_cons <- sapply(unique(cosmic_vep$Uploaded_variation), cons_fea, input_data = cosmic_vep, vep_levels = vep_levels)
cosmic_vep_cons <- data.frame(t(cosmic_vep_cons))

genome_vep_cons <- sapply(unique(genome_vep$Uploaded_variation), cons_fea, input_data = genome_vep, vep_levels = vep_levels)
genome_vep_cons <- data.frame(t(genome_vep_cons))

save(cosmic_vep_cons, genome_vep_cons, file = "dataset_v1.3_vepCons.RData")

gencode_chrx <- read.table("gencode.v31.annotation.chrX.gtf", sep = "\t", stringsAsFactors = F)
input_info <- gencode_chrx$V9[[55]]

splitINFO <- function(input_info) {
        INFO <- strsplit(input_info, split = "; ")[[1]]
        idx_tName <- which(gregexpr("transcript_type ", INFO, fixed = T) != -1)
        idx_gName <- which(gregexpr("gene_type ", INFO, fixed = T) != -1)
        
        if(length(idx_gName) != 0) {
                gene_type <- INFO[idx_gName]
                gene_type <- gsub("gene_type ", "", gene_type)
        } else gene_type <- "-"
        
        if(length(idx_tName) != 0) {
                transcript_type <- INFO[idx_tName]
                transcript_type <- gsub("transcript_type ", "", transcript_type)
        } else transcript_type <- "-"
        
        output <- c(gene_type = gene_type, transcript_type = transcript_type)
        output
}

type_INFO <- sapply(gencode_chrx$V9, splitINFO, USE.NAMES = F)

type_INFO_df <- data.frame(t(type_INFO), stringsAsFactors = F)

gencode_chrx <- cbind(gencode_chrx, type_INFO_df)

input_POS <- GenomeData_hg38_window$POS[[10]]
input_POS <- 1288808

annote_region <- function(input_POS, gencode_chrx = gencode_chrx) {
        new_df <- gencode_chrx[((input_POS >= gencode_chrx$V4) & (input_POS <= gencode_chrx$V5)),]
}


cd_cons <- rbind(cosmicCD_consv_vep_sampled[,14:49], genomeCD_consv_vep[,13:48])
nc_cons <- rbind(cosmicNC_consv_vep[,14:49], genomeNC_consv_vep_sampled[,13:48])

cd_cons <- cosmic_chrX_all_coding_VEP[,2:37]
nc_cons <- GenomeProj_VEP_nc[,2:37]

tmp <- sort(apply(cd_cons, 2, sum), decreasing = T)
tmp_df <- data.frame(count = tmp, percentage = round(tmp / 203834 * 100, digits = 2), 
                     row.names = sapply(names(tmp), gsub, pattern = "_", replacement = " ", USE.NAMES = F))
write.csv(tmp_df, file = "cd_consequence_type.csv", quote = F)

###### distribution ######

cosmic_vep <- read.table("vep_CosmicNonCodingVariants_ChrX", sep = "\t", header = T, stringsAsFactors = F)

all_len <- length(unique(cosmic_vep$Uploaded_variation))

# idx <- seq(1, 203834, 6795)
idx <- seq(1, all_len, ceiling(all_len / 30))

unique_names <- unique(cosmic_vep$Uploaded_variation)
res <- lapply(idx, function(x) {
        bound <- ifelse(x + (ceiling(all_len / 30) - 1) < all_len, (x + ceiling(all_len / 30) - 1), all_len)
        res <- unique_names[x:bound]
})

# length(unique(unlist(res)))

AA_flag <- function(aminoAcidINFO, AA.LETTERS = sort(c("G", "V", "L", "F", "P", "I", "A", "M", 
                                                       "T", "S", "Y", "N", "Q", "W", "H", "K", 
                                                       "R", "E", "D", "C"))) {
        aa_vec <- strsplit(aminoAcidINFO, split = "/")[[1]]
        
        if (length(aa_vec) == 2) {
                if (aa_vec[[1]] == "*") REF_vec <- rep(1, 20) else REF_vec <- as.numeric(AA.LETTERS %in% aa_vec[[1]])
                if (aa_vec[[2]] == "*") ALT_vec <- rep(1, 20) else ALT_vec <- as.numeric(AA.LETTERS %in% aa_vec[[2]])
        } else if (aminoAcidINFO %in% AA.LETTERS) {
                if (aminoAcidINFO == "*") REF_vec <- rep(1, 20) else REF_vec <- as.numeric(AA.LETTERS %in% aminoAcidINFO[[1]])
                ALT_vec <- REF_vec
        } else {
                REF_vec <- rep(0, 20)
                ALT_vec <- REF_vec
        }
        out_vec <- c(REF_vec, ALT_vec)
        out_vec
}

all_flag_name <- function(idx_num, idx_name, input_data, vep_levels = vep_levels, 
                          AA.LETTERS = sort(c("G", "V", "L", "F", "P", "I", "A", "M", 
                                              "T", "S", "Y", "N", "Q", "W", "H", "K", 
                                              "R", "E", "D", "C"))) {
        idx <- idx_name[idx_num]
        if(idx_num %% 1000 == 0) message(idx_num, " / ", length(idx_name), " ",Sys.time())
        
        new_data <- input_data[input_data$Uploaded_variation == idx,]
        # res <- sapply(new_data$Amino_acids, AA_flag, AA.LETTERS = AA.LETTERS)
        # res <- apply(res, 1, max)
        # # names(res) <- c(paste0("REF_", AA.LETTERS), paste0("ALT_", AA.LETTERS))
        # aa_res <- res
        
        res <- lapply(new_data$Consequence, function(x) {
                out <- strsplit(x, ",")[[1]]
        })
        # browser()
        res <- unique(unlist(res))
        res <- table(factor(res, levels = vep_levels))
        cons_res <- c(num = nrow(new_data), res)
        
        # output <- c(cons_res, aa_res)
}

wrap_func <- function(idx_name, input_data = cosmic_vep) {
        
        vep_levels <- strsplit(("synonymous_variant, non_coding_transcript_exon_variant, missense_variant, stop_gained, frameshift_variant, incomplete_terminal_codon_variant, coding_sequence_variant, downstream_gene_variant, splice_region_variant, intron_variant, stop_retained_variant, upstream_gene_variant, splice_donor_variant, non_coding_transcript_variant, splice_acceptor_variant, 3_prime_UTR_variant, NMD_transcript_variant, 5_prime_UTR_variant, start_lost, inframe_deletion, stop_lost, inframe_insertion, protein_altering_variant, intergenic_variant, mature_miRNA_variant, start_retained_variant, transcript_ablation, transcript_amplification, TFBS_ablation, TFBS_amplification, TF_binding_site_variant, regulatory_region_ablation, regulatory_region_amplification, feature_elongation, regulatory_region_variant, feature_truncation"), 
                               ", ")[[1]]
        
        cosmic_vep_all <- sapply(1:length(idx_name), 
                                 all_flag_name, idx_name = idx_name, 
                                 input_data = input_data, vep_levels = vep_levels,
                                 AA.LETTERS = sort(c("G", "V", "L", "F", "P", "I", "A", "M", 
                                                     "T", "S", "Y", "N", "Q", "W", "H", "K", 
                                                     "R", "E", "D", "C")))
        cosmic_vep_all
}

cl <- parallel::makeCluster(30, outfile = "")
parallel::clusterExport(cl, varlist = c("AA_flag", "all_flag_name"))

# cosmic_vep <- read.table("vep_CosmicNonCodingVariants_ChrX", sep = "\t", header = T, stringsAsFactors = F)

res_noncoding <- parallel::parSapply(cl, res, wrap_func, input_data = cosmic_vep)
parallel::stopCluster(cl)
save(res_noncoding, file = "res_noncoding.RData")
res_noncoding_df <- do.call("cbind", res_noncoding)
res_noncoding_df <- data.frame(t(res_noncoding_df))
save(res_noncoding_df, file = "res_noncoding_df.RData")

####

cosmic_vep <- read.table("vep_GenomeProj_chrX_all.vcf", sep = "\t", header = T, stringsAsFactors = F)

all_len <- length(unique(cosmic_vep$Uploaded_variation))

# idx <- seq(1, 203834, 6795)
idx <- seq(1, all_len, ceiling(all_len / 30))

unique_names <- unique(cosmic_vep$Uploaded_variation)
res <- lapply(idx, function(x) {
        bound <- ifelse(x + (ceiling(all_len / 30) - 1) < all_len, (x + ceiling(all_len / 30) - 1), all_len)
        res <- unique_names[x:bound]
})

# length(unique(unlist(res))) == all_len

AA_flag <- function(aminoAcidINFO, AA.LETTERS = sort(c("G", "V", "L", "F", "P", "I", "A", "M", 
                                                       "T", "S", "Y", "N", "Q", "W", "H", "K", 
                                                       "R", "E", "D", "C"))) {
        aa_vec <- strsplit(aminoAcidINFO, split = "/")[[1]]
        
        if (length(aa_vec) == 2) {
                if (aa_vec[[1]] == "*") REF_vec <- rep(1, 20) else REF_vec <- as.numeric(AA.LETTERS %in% aa_vec[[1]])
                if (aa_vec[[2]] == "*") ALT_vec <- rep(1, 20) else ALT_vec <- as.numeric(AA.LETTERS %in% aa_vec[[2]])
        } else if (aminoAcidINFO %in% AA.LETTERS) {
                if (aminoAcidINFO == "*") REF_vec <- rep(1, 20) else REF_vec <- as.numeric(AA.LETTERS %in% aminoAcidINFO[[1]])
                ALT_vec <- REF_vec
        } else {
                REF_vec <- rep(0, 20)
                ALT_vec <- REF_vec
        }
        out_vec <- c(REF_vec, ALT_vec)
        out_vec
}

all_flag_name <- function(idx_num, idx_name, input_data, vep_levels = vep_levels, 
                          AA.LETTERS = sort(c("G", "V", "L", "F", "P", "I", "A", "M", 
                                              "T", "S", "Y", "N", "Q", "W", "H", "K", 
                                              "R", "E", "D", "C"))) {
        idx <- idx_name[idx_num]
        if(idx_num %% 1000 == 0) message(idx_num, " / ", length(idx_name), " ",Sys.time())
        
        new_data <- input_data[input_data$Uploaded_variation == idx,]
        res <- sapply(new_data$Amino_acids, AA_flag, AA.LETTERS = AA.LETTERS)
        res <- apply(res, 1, max)
        names(res) <- c(paste0("REF_", AA.LETTERS), paste0("ALT_", AA.LETTERS))
        aa_res <- res
        
        res <- lapply(new_data$Consequence, function(x) {
                out <- strsplit(x, ",")[[1]]
        })
        # browser()
        res <- unique(unlist(res))
        res <- table(factor(res, levels = vep_levels))
        cons_res <- c(num = nrow(new_data), res)
        
        output <- c(cons_res, aa_res)
}

wrap_func <- function(idx_name, input_data = cosmic_vep) {
        
        vep_levels <- strsplit(("synonymous_variant, non_coding_transcript_exon_variant, missense_variant, stop_gained, frameshift_variant, incomplete_terminal_codon_variant, coding_sequence_variant, downstream_gene_variant, splice_region_variant, intron_variant, stop_retained_variant, upstream_gene_variant, splice_donor_variant, non_coding_transcript_variant, splice_acceptor_variant, 3_prime_UTR_variant, NMD_transcript_variant, 5_prime_UTR_variant, start_lost, inframe_deletion, stop_lost, inframe_insertion, protein_altering_variant, intergenic_variant, mature_miRNA_variant, start_retained_variant, transcript_ablation, transcript_amplification, TFBS_ablation, TFBS_amplification, TF_binding_site_variant, regulatory_region_ablation, regulatory_region_amplification, feature_elongation, regulatory_region_variant, feature_truncation"), 
                               ", ")[[1]]
        
        cosmic_vep_all <- sapply(1:length(idx_name), 
                                 all_flag_name, idx_name = idx_name, 
                                 input_data = input_data, vep_levels = vep_levels,
                                 AA.LETTERS = sort(c("G", "V", "L", "F", "P", "I", "A", "M", 
                                                     "T", "S", "Y", "N", "Q", "W", "H", "K", 
                                                     "R", "E", "D", "C")))
        cosmic_vep_all
}

cl <- parallel::makeCluster(30, outfile = "")
parallel::clusterExport(cl, varlist = c("AA_flag", "all_flag_name"))

# cosmic_vep <- read.table("vep_CosmicNonCodingVariants_ChrX", sep = "\t", header = T, stringsAsFactors = F)

GenomeProj_VEP_Cons_AA <- parallel::parSapply(cl, res, wrap_func, input_data = cosmic_vep)
parallel::stopCluster(cl)
save(GenomeProj_VEP_Cons_AA, file = "GenomeProj_VEP_Cons_AA.RData")
GenomeProj_VEP_Cons_AA_df <- do.call("cbind", GenomeProj_VEP_Cons_AA)
GenomeProj_VEP_Cons_AA_df <- data.frame(t(GenomeProj_VEP_Cons_AA_df))
save(GenomeProj_VEP_Cons_AA_df, file = "GenomeProj_VEP_Cons_AA_df.RData")

#########

chk_aa_sum <- apply(GenomeProj_VEP_Cons_AA_df[,38:77], 1, sum)

GenomeProj_VEP_nc <- GenomeProj_VEP_Cons_AA_df[chk_aa_sum == 0,]
GenomeProj_VEP_cd <- GenomeProj_VEP_Cons_AA_df[chk_aa_sum != 0,]
save(GenomeProj_VEP_cd, GenomeProj_VEP_nc, file = "GenomeProj_VEP_cd_nc.RData")

rm(GenomeProj_VEP_Cons_AA_df, chk_aa_sum)

cosmic_num_cd <- cosmic_chrX_all_coding_VEP$num
cosmic_num_nc <- cosmic_chrX_all_noncoding_VEP$num
genome_num_cd <- GenomeProj_VEP_cd$num
genome_num_nc <- GenomeProj_VEP_nc$num

cosmic_cd <- data.frame(Num = cosmic_num_cd, Source = "COSMIC", Region = "Coding")
cosmic_nc <- data.frame(Num = cosmic_num_nc, Source = "COSMIC", Region = "NonCoding")

genome_cd <- data.frame(Num = genome_num_cd, Source = "1000 Genomes", Region = "Coding")
genome_nc <- data.frame(Num = genome_num_nc, Source = "1000 Genomes", Region = "NonCoding")

num_val <- rbind(cosmic_cd, cosmic_nc, genome_cd, genome_nc)

library(ggplot2)

ggplot(data = num_val, aes(x = Region, y = Num)) + 
        geom_boxplot(alpha = I(0.6), aes(fill = Source)) +
        scale_y_sqrt(limits = c(1, 60), breaks = c(2,4,6,8,10,15,25,40,60)) +
        coord_flip()

cd_cons <- read.csv("summary_cd_vep_consequence.csv")

cd_cons$Variant_type <- factor(cd_cons$Variant_type, 
                               levels = unique(as.character(cd_cons$Variant_type[order(cd_cons$Percentage, decreasing = T)])))

ggplot(data = cd_cons, aes(x = Variant_type, y = Percentage, fill = Source)) +
        geom_bar(stat = "identity", position = position_dodge(), colour = "black") +
        scale_y_continuous(breaks = seq(0, 100, 10)) +
        theme(axis.text.x = element_text(angle = 10, vjust = 1, hjust = 1, size = 14),
              axis.title.x = element_text(size = 14),
              axis.title.y = element_text(size = 14),
              axis.text.y  = element_text(size = 14),
              legend.title = element_text(size = 14),
              legend.text  = element_text(size = 14)) # 16 * 5

nc_cons <- read.csv("summary_nc_vep_consequence.csv")

nc_cons$Variant_type <- factor(nc_cons$Variant_type, 
                               levels = unique(as.character(nc_cons$Variant_type[order(nc_cons$Percentage, decreasing = T)])))

ggplot(data = nc_cons, aes(x = Variant_type, y = Percentage, fill = Source)) +
        geom_bar(stat = "identity", position = position_dodge(), colour = "black") +
        scale_y_continuous(breaks = seq(0, 100, 10)) +
        theme(axis.text.x = element_text(angle = 10, vjust = 1, hjust = 1, size = 14),
              axis.title.x = element_text(size = 14),
              axis.title.y = element_text(size = 14),
              axis.text.y  = element_text(size = 14),
              legend.title = element_text(size = 14),
              legend.text  = element_text(size = 14)) # 16 * 5


data_nc <- data.frame(rbind(cbind(data = "COSMIC_cd", num = cosmic_chrX_all_coding_VEP$num),
                            cbind(data = "1000G_cd", num = GenomeProj_VEP_cd$num),
                            cbind(data = "COSMIC_nc", num = cosmic_chrX_all_noncoding_VEP$num),
                            cbind(data = "1000G_nc", num = GenomeProj_VEP_nc$num)), stringsAsFactors = F)
data_nc$num <- as.numeric(data_nc$num)
library(ggplot2)

ggplot(data_nc, aes(x = num, color = data, fill = data)) + 
        geom_density(bw = 1, alpha = 0.1, size = 1)  +
        # coord_cartesian(ylim = c(0, 10)) + 
        scale_x_continuous(breaks = seq(0, 90, 5)) +
        theme(axis.text.x = element_text(size = 14),
              axis.title.x = element_text(size = 14),
              axis.title.y = element_text(size = 14),
              axis.text.y  = element_text(size = 14),
              legend.title = element_text(size = 14),
              legend.text  = element_text(size = 14)) # 3 * 12 in

#######

set.seed(1)
cd_idx <- sample(432, 345)
nc_idx <- sample(5260, 4000)

train_NC <- rbind(cosmicNC_consv_vep[(1:5260 %in% nc_idx),-c(1:6)], 
                  genomeNC_consv_vep[(1:5260 %in% nc_idx),-c(1:6)])
test_NC  <- rbind(cosmicNC_consv_vep[(!1:5260 %in% nc_idx),-c(1:6)], 
                  genomeNC_consv_vep[(!1:5260 %in% nc_idx),-c(1:6)])

train_CD <- rbind(cosmicCD_consv_vep[(1:432 %in% cd_idx),-c(1:6)], 
                  genomeCD_consv_vep[(1:432 %in% cd_idx),-c(1:6)])
test_CD  <- rbind(cosmicCD_consv_vep[(!1:432 %in% cd_idx),-c(1:6)], 
                  genomeCD_consv_vep[(!1:432 %in% cd_idx),-c(1:6)])
save.image(file = "dataset_v1.6.1_consv_vep_sampled_trainTest.RData")

mod_cd_consv <- randomForest_tune(datasets = list(train_CD[,1:5]),
                                  ntree.range = c(50,100,200,300,400,500,600,800,1000),
                                  return.model = T, parallel.cores = -1, seed = 1)
mod_cd_vep_cons <- randomForest_tune(datasets = list(train_CD[,c(1, 6:42)]),
                                     ntree.range = c(50,100,200,300,400,500,600,800,1000),
                                     return.model = T, parallel.cores = -1, seed = 1)
mod_cd_vep_aa <- randomForest_tune(datasets = list(train_CD[,c(1, 43:82)]),
                                   ntree.range = c(50,100,200,300,400,500,600,800,1000),
                                   return.model = T, parallel.cores = -1, seed = 1)
mod_cd_vep_all <- randomForest_tune(datasets = list(train_CD[,c(1, 6:82)]),
                                    ntree.range = c(50,100,200,300,400,500,600,800,1000),
                                    return.model = T, parallel.cores = -1, seed = 1)
mod_cd_consv_vepCons <- randomForest_tune(datasets = list(train_CD[,1:42]),
                                          ntree.range = c(50,100,200,300,400,500,600,800,1000),
                                          return.model = T, parallel.cores = -1, seed = 1)
mod_cd_consv_vepAll <- randomForest_tune(datasets = list(train_CD),
                                         ntree.range = c(50,100,200,300,400,500,600,800,1000),
                                         return.model = T, parallel.cores = -1, seed = 1)

mod_nc_consv <- randomForest_tune(datasets = list(train_NC[,1:5]),
                                  ntree.range = c(50,100,200,300,400,500,600,800,1000),
                                  return.model = T, parallel.cores = -1, seed = 1)
mod_nc_vep_cons <- randomForest_tune(datasets = list(train_NC[,c(1, 6:42)]),
                                     ntree.range = c(50,100,200,300,400,500,600,800,1000),
                                     return.model = T, parallel.cores = -1, seed = 1)
mod_nc_consv_vepCons <- randomForest_tune(datasets = list(train_NC[,1:42]),
                                          ntree.range = c(50,100,200,300,400,500,600,800,1000),
                                          return.model = T, parallel.cores = -1, seed = 1)

library(randomForest)

res_cd_consv         <- predict(mod_cd_consv, test_CD, type = "prob")
res_cd_vep_cons      <- predict(mod_cd_vep_cons, test_CD, type = "prob")
res_cd_vep_aa        <- predict(mod_cd_vep_aa, test_CD, type = "prob")
res_cd_vep_all       <- predict(mod_cd_vep_all, test_CD, type = "prob")
res_cd_consv_vepCons <- predict(mod_cd_consv_vepCons, test_CD, type = "prob")
res_cd_consv_vepAll  <- predict(mod_cd_consv_vepAll, test_CD, type = "prob")

res_nc_consv         <- predict(mod_nc_consv, test_NC, type = "prob")
res_nc_vep_cons      <- predict(mod_nc_vep_cons, test_NC, type = "prob")
res_nc_consv_vepCons <- predict(mod_nc_consv_vepCons, test_NC, type = "prob")

for (i in ls(pattern = "res_nc")) {
        res <- get(i)
        res_data <- pROC::roc(test_NC$label,  res[,1])
        assign(paste0("roc_", gsub("res_", "", i)), res_data)
        print(paste0(i, ": ", round(pROC::auc(test_NC$label,  res[,1]), 4)))
}

library(pROC)
pROC::plot.roc(cpc.roc,  print.auc = F, legacy.axes = T,
               print.thres = F, col = "dodgerblue4")
pROC::plot.roc(cpat.roc, print.auc = F, add = T, legacy.axes = T,
               print.thres = F, col = "cyan3")
pROC::plot.roc(cnci.roc, print.auc = F, add = T, legacy.axes = T,
               print.thres = F, col = "maroon2")
pROC::plot.roc(plek.roc, print.auc = F, add = T, legacy.axes = T,
               print.thres = F, col = "green3")
pROC::plot.roc(cpc2.roc, print.auc = F, add = T, legacy.axes = T,
               print.thres = F, col = "plum3")
pROC::plot.roc(tool.roc, print.auc = F, add = T, legacy.axes = T,
               print.thres = F, col = "firebrick2")

pROC::ggroc(data = list(consv = roc_cd_consv, consv_vepAll = roc_cd_consv_vepAll, 
                        consv_vepCons = roc_cd_consv_vepCons, vep_AA = roc_cd_vep_aa, 
                        vep_all = roc_cd_vep_all, vep_Cons = roc_cd_vep_cons), 
            legacy.axes = T, size = 1) 

pROC::ggroc(data = list(consv = roc_nc_consv,
                        consv_vepCons = roc_nc_consv_vepCons,
                        vep_Cons = roc_nc_vep_cons), 
            legacy.axes = T, size = 1)  # 6 * 8

######

annote_type <- function(input_type, input_POS, input_df, threshold) {
        new_df <- input_df[input_df$V3 == input_type,]
        distance_range_start <- abs(input_POS - new_df$V4)
        distance_range_end   <- abs(input_POS - new_df$V5)
        distance_min <- min(c(distance_range_start, distance_range_end))
        if (distance_min <= threshold) {
                distance <- 1 / (distance_min + 1)
        } else distance <- 0
        distance
}

annote_region <- function(input_POS, gtf_data, threshold) {
        distance_all <- c(CDS = 0, exon = 0, gene = 0, start_codon = 0, stop_codon = 0, transcript = 0, UTR = 0)
        new_df_all <- gtf_data[((input_POS >= gtf_data$V4) & (input_POS <= gtf_data$V5)),]
        if (nrow(new_df_all) != 0) {
                res <- sapply(unique(new_df_all$V3), annote_type, input_POS = input_POS, 
                              input_df = new_df_all, threshold = threshold)
                for (i in 1:length(res)) {
                        distance_all[names(distance_all) == names(res)[i]] <- res[i]
                }
        }
        distance_all
}

annote_NC <- function(all_input_POS, gtf_data = gencode_chrx, threshold, label = NULL, parallel.cores = 2, cl = NULL) {
        if (is.null(cl)) {
                cl <- parallel::makeCluster(parallel.cores)
                parallel::clusterExport(cl, "annote_type")
        }
        
        distance_test <- parallel::parSapply(cl, all_input_POS, annote_region, 
                                             gtf_data = gtf_data, threshold = threshold)
        
        if (is.null(cl)) {
                parallel::stopCluster(cl)
        }
        
        distance_test_df <- data.frame(t(distance_test), stringsAsFactors = F)
        if (!is.null(label)) distance_test_df <- cbind(label = label, distance_test_df, stringsAsFactors = F)
        distance_test_df
}

# distanceDataset_10000 <- annote_NC(all_input_POS = c(cosmicNC_consv_vep$POS, genomeNC_consv_vep$POS), 
#                                    gtf_data = gencode_chrx, threshold = 10000, parallel.cores = 2, 
#                                    label = c(cosmicNC_consv_vep$label, genomeNC_consv_vep$label))

cl <- parallel::makeCluster(10)
parallel::clusterExport(cl, "annote_type")

for (threshold_input in c(1, 10, 50, 100, 200, 300, 500, 700, 1000, 1500, 2000, 2500, 
                         3000, 3500, seq(4000, 10000, 1000))) {
        message("Proccessing: ", threshold_input, "     ", Sys.time())
        
        distanceDataset <- annote_NC(all_input_POS = c(cosmicCD_consv_vep$POS, genomeCD_consv_vep$POS), 
                                     gtf_data = gencode_chrx, threshold = threshold_input, cl = cl, 
                                     label = c(as.character(cosmicCD_consv_vep$label), 
                                               as.character(genomeCD_consv_vep$label)))
        distanceDataset$label <- as.factor(distanceDataset$label)
        assign(paste0("distanceFeatureSet_", threshold_input), distanceDataset)
        # save(list = ls(pattern = "distanceFeatureSet_"), file = "distanceFeatureSet_cd.RData")
}
parallel::stopCluster(cl)

# save.image(file = "tmp_distance.RData", compress = "xz")

# for (data_name in ls(pattern = "distanceFeatureSet_")) {
#         dataset_input <- get(data_name)
#         dataset_input$label <- as.factor(dataset_input$label)
#         assign(data_name, dataset_input)
# }

for (data_name in ls(pattern = "distanceFeatureSet_")) {
        message("Proccessing: ", data_name, "     ", Sys.time())
        dataset_input <- get(data_name)
        res_distance <- randomForest_tune(datasets = list(dataset_input), label.col = 1,
                                          positive.class = "Positive", folds.num = 10,
                                          ntree.range = c(50, 100, 200, 300, 500, 700, 1000, 1500),
                                          seed = 1, return.model = F, parallel.cores = 10)
        assign(paste0("res_", data_name), res_distance)
        save(list = ls(pattern = "res_distanceFeatureSet_"), file = "res_distance_cd.RData")
}

# res_distance <- randomForest_tune(datasets = list(distanceFeatureSet_1), label.col = 1,
#                                   positive.class = "Positive", folds.num = 10,
#                                   ntree.range = c(50, 100, 200, 300, 500, 700, 1000, 1500),
#                                   seed = 1, return.model = F, parallel.cores = 10)

perf_distance <- c()
for (i in ls(pattern = "res_distance")) {
        res_data <- get(i)
        perf_out <- res_data$performance[row.names(res_data$performance) == paste0("ntree_", res_data$ntree),]
        perf_out$WindowSize <- as.numeric(gsub("res_distanceFeatureSet_", "", i))
        perf_distance <- rbind (perf_distance, perf_out)
}

perf_distance$Region = "Coding"

row.names(perf_distance) <- NULL
perf_distance_cd <- perf_distance
save(perf_distance_cd, file = "perf_distance_cd.RData")

######

load("tmp_distance.RData")

cl <- parallel::makeCluster(20)
parallel::clusterExport(cl, "annote_type")

for (threshold_input in seq(50000, 100000, 5000)) {
        message("Proccessing: ", threshold_input, "     ", Sys.time())
        
        distanceDatasetCD <- annote_NC(all_input_POS = c(cosmicCD_consv_vep$POS, genomeCD_consv_vep$POS), 
                                       gtf_data = gencode_chrx, threshold = threshold_input, cl = cl, 
                                       label = c(as.character(cosmicCD_consv_vep$label), 
                                                 as.character(genomeCD_consv_vep$label)))
        distanceDatasetNC <- annote_NC(all_input_POS = c(cosmicNC_consv_vep$POS, genomeNC_consv_vep$POS), 
                                       gtf_data = gencode_chrx, threshold = threshold_input, cl = cl, 
                                       label = c(as.character(cosmicNC_consv_vep$label), 
                                                 as.character(genomeNC_consv_vep$label)))
        
        distanceDatasetCD$label <- as.factor(distanceDatasetCD$label)
        distanceDatasetNC$label <- as.factor(distanceDatasetNC$label)
        
        assign(paste0("distanceFeatureSet_cd_", threshold_input), distanceDatasetCD)
        assign(paste0("distanceFeatureSet_nc_", threshold_input), distanceDatasetNC)
        save(list = ls(pattern = "distanceFeatureSet_cd_"), file = "distanceFeatureSet_cd_part4.RData")
        save(list = ls(pattern = "distanceFeatureSet_nc_"), file = "distanceFeatureSet_nc_part4.RData")
}
parallel::stopCluster(cl)

for (data_name in ls(pattern = "distanceFeatureSet_cd_")) {
        message("Proccessing: ", data_name, "     ", Sys.time())
        dataset_input <- get(data_name)
        res_distance <- randomForest_tune(datasets = list(dataset_input), label.col = 1,
                                          positive.class = "Positive", folds.num = 10,
                                          ntree.range = c(50, 100, 200, 300, 500, 700, 1000, 1500),
                                          seed = 1, return.model = F, parallel.cores = 10)
        assign(paste0("res_", data_name), res_distance)
        save(list = ls(pattern = "res_distanceFeatureSet_"), file = "res_distance_cd_part4.RData")
}

for (data_name in ls(pattern = "distanceFeatureSet_nc_")) {
        message("Proccessing: ", data_name, "     ", Sys.time())
        dataset_input <- get(data_name)
        res_distance <- randomForest_tune(datasets = list(dataset_input), label.col = 1,
                                          positive.class = "Positive", folds.num = 10,
                                          ntree.range = c(50, 100, 200, 300, 500, 700, 1000, 1500),
                                          seed = 1, return.model = F, parallel.cores = 10)
        assign(paste0("res_", data_name), res_distance)
        save(list = ls(pattern = "res_distanceFeatureSet_"), file = "res_distance_nc_part4.RData")
}

#

format_res <- function(input_pattern, Region = c("Coding", "NonCoding")) {
        perf_distance <- c()
        file_list <- ls(pattern = input_pattern, envir = globalenv())
        print(file_list)
        for (i in file_list) {
                res_data <- get(i)
                perf_out <- res_data$performance[row.names(res_data$performance) == paste0("ntree_", res_data$ntree),]
                perf_out$WindowSize <- as.numeric(gsub(input_pattern, "", i))
                perf_distance <- rbind(perf_distance, perf_out)
        }
        
        perf_distance$Region = Region
        
        row.names(perf_distance) <- NULL
        perf_distance
}

res_distance_cd_part4 <- format_res(input_pattern = "res_distanceFeatureSet_cd_", Region = "Coding")
res_distance_nc_part4 <- format_res(input_pattern = "res_distanceFeatureSet_nc_", Region = "NonCoding")

save(res_distance_cd_part4, res_distance_nc_part4, file = "perf_distance_cd_nc_part4.RData")

res_distance_cd <- rbind(perf_distance_cd, res_distance_cd_part2, res_distance_cd_part3, res_distance_cd_part4)
res_distance_nc <- rbind(perf_distance_nc, res_distance_nc_part2, res_distance_nc_part3, res_distance_nc_part4)

res_distance_cd <- res_distance_cd[order(res_distance_cd$WindowSize),]
res_distance_nc <- res_distance_nc[order(res_distance_nc$WindowSize),]

write.csv(res_distance_cd, file = "res_distance_cd.csv", row.names = F)
write.csv(res_distance_nc, file = "res_distance_nc.csv", row.names = F)
save(res_distance_cd, res_distance_nc, file = "res_distance_nc_cd.RData")

library(ggplot2)

perf_distance <- read.table("res_distance_all_ggplot.csv", sep = ",", header = T)

perf_distance <- perf_distance[perf_distance$Metric %in% c("Accuracy", "F-Measure"),]
perf_distance$Window_Size <- as.factor(perf_distance$Window_Size)

ggplot(data = perf_distance, aes(x = Window_Size, y = Performance, group = Group, 
                                 colour = Region)) +
        geom_line(size = 1, aes(linetype = Metric)) + 
        scale_y_continuous(breaks = seq(0, 1, 0.1)) +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8),
              axis.title.x = element_text(size = 12),
              axis.title.y = element_text(size = 12),
              axis.text.y  = element_text(size = 12),
              legend.title = element_text(size = 12),
              legend.text  = element_text(size = 12)) # 3 * 12 in


###


annote_type <- function(input_type, input_POS, input_df) {
        new_df <- input_df[input_df$V3 == input_type,]
        distance_range_start <- abs(input_POS - new_df$V4)
        distance_range_end   <- abs(input_POS - new_df$V5)
        distance_min <- min(c(distance_range_start, distance_range_end))
        distance_min
}

annote_region <- function(input_POS, gtf_data) {
        distance_all <- c(CDS = 0, exon = 0, gene = 0, start_codon = 0, stop_codon = 0, transcript = 0, UTR = 0)
        new_df_all <- gtf_data[((input_POS >= gtf_data$V4) & (input_POS <= gtf_data$V5)),]
        if (nrow(new_df_all) != 0) {
                res <- sapply(unique(new_df_all$V3), annote_type, input_POS = input_POS, 
                              input_df = new_df_all)
                for (i in 1:length(res)) {
                        distance_all[names(distance_all) == names(res)[i]] <- res[i]
                }
        }
        distance_all
}

annote_NC <- function(all_input_POS, gtf_data = gencode_chrx, label = NULL, parallel.cores = 2, cl = NULL) {
        if (is.null(cl)) {
                cl <- parallel::makeCluster(parallel.cores)
                parallel::clusterExport(cl, "annote_type")
        }
        
        distance_test <- parallel::parSapply(cl, all_input_POS, annote_region, 
                                             gtf_data = gtf_data)
        
        if (is.null(cl)) {
                parallel::stopCluster(cl)
        }
        
        distance_test_df <- data.frame(t(distance_test), stringsAsFactors = F)
        if (!is.null(label)) distance_test_df <- cbind(label = label, distance_test_df, stringsAsFactors = F)
        distance_test_df
}


distanceDatasetCD <- annote_NC(all_input_POS = c(cosmicCD_consv_vep$POS, genomeCD_consv_vep$POS), 
                               gtf_data = gencode_chrx , cl = cl, 
                               label = c(as.character(cosmicCD_consv_vep$label), 
                                         as.character(genomeCD_consv_vep$label)))
distanceDatasetNC <- annote_NC(all_input_POS = c(cosmicNC_consv_vep$POS, genomeNC_consv_vep$POS), 
                               gtf_data = gencode_chrx, cl = cl, 
                               label = c(as.character(cosmicNC_consv_vep$label), 
                                         as.character(genomeNC_consv_vep$label)))


save(distanceDatasetNC, distanceDatasetCD, file = "distance_distribution.RData")

distance_range <- rbind(distanceDatasetCD, distanceDatasetNC)
rm(distanceDatasetNC, distanceDatasetCD)

quantile(distance_range_NC, probs = seq(0.1, 1, 0.1))
quantile(distance_range_CD, probs = seq(0.1, 1, 0.1))

tmp <- sapply(distance_range, function(x) {
        x[x == 0] <- NA
        x
})

tmp <- data.frame(tmp, stringsAsFactors = F)

write.csv(tmp, file = "tmp.csv", quote = F, row.names = F)

tmp <- read.csv("tmp.csv")

library(ggplot2)

ggplot(data = tmp, aes(x = GenomeFeature, y = Distance)) + 
        geom_boxplot(alpha = I(0.6), aes(fill = Region)) +
        scale_y_continuous(limits = c(0, 25000))

########

perf_cv_consv_dist_cd <- randomForest_tune(datasets = list(data_cd[,6:17]),
                                           ntree.range = c(50,100,200,300,400,500,800,1000,1500),
                                           return.model = F, parallel.cores = -1, seed = 1)
perf_cv_consv_dist_nc <- randomForest_tune(datasets = list(data_nc[,6:17]),
                                           ntree.range = c(50,100,200,300,400,500,800,1000,1500),
                                           return.model = F, parallel.cores = -1, seed = 1)
save(perf_cv_consv_dist_cd, perf_cv_consv_dist_nc, file = "perf_cv_consv_dist_cd_nc.RData")

write.csv(perf_cv_consv_dist_cd$performance, file = "perf_cv_consv_dist_cd.csv", quote = F)
write.csv(perf_cv_consv_dist_nc$performance, file = "perf_cv_consv_dist_nc.csv", quote = F)

perf_cv_consv_dist_cd <- randomForest_tune(datasets = list(data_cd[,6:17]),
                                           ntree.range = c(50,100,200,300,400,500,800,1000,1500),
                                           return.model = F, parallel.cores = -1, seed = 1)

perf_cv_consv_dist_vepCons_nc <- randomForest_tune(datasets = list(data_nc[,6:54]),
                                                   ntree.range = c(50,100,200,300,400,500,800,1000,1500),
                                                   return.model = F, parallel.cores = -1, seed = 1)
write.csv(perf_cv_consv_dist_vepCons_nc$performance, file = "perf_cv_consv_dist_vepCons_nc.csv", quote = F)

perf_cv_consv_dist_vepCons_cd <- randomForest_tune(datasets = list(data_cd[,6:54]),
                                                   ntree.range = c(50,100,200,300,400,500,800,1000,1500),
                                                   return.model = F, parallel.cores = -1, seed = 1)
write.csv(perf_cv_consv_dist_vepCons_cd$performance, file = "perf_cv_consv_dist_vepCons_cd.csv", quote = F)

pos_nc <- data_nc[data_nc$label == "Positive",]
neg_nc <- data_nc[data_nc$label != "Positive",]

pos_cd <- data_cd[data_cd$label == "Positive",]
neg_cd <- data_cd[data_cd$label != "Positive",]

train_NC <- rbind(pos_nc[(1:5260 %in% nc_idx),-c(1:5)], 
                  neg_nc[(1:5260 %in% nc_idx),-c(1:5)])

test_NC  <- rbind(pos_nc[(!1:5260 %in% nc_idx),-c(1:5)], 
                  neg_nc[(!1:5260 %in% nc_idx),-c(1:5)])

train_CD <- rbind(pos_cd[(1:432 %in% cd_idx),-c(1:5)], 
                  neg_cd[(1:432 %in% cd_idx),-c(1:5)])
test_CD  <- rbind(pos_cd[(!1:432 %in% cd_idx),-c(1:5)], 
                  neg_cd[(!1:432 %in% cd_idx),-c(1:5)])

save.image(file = "dataset_v1.7.1_consv_dist_vep_sampled_trainTest.RData")

load("dataset_v1.7.1_consv_dist_vep_sampled_trainTest.RData")

mod_cd_dist <- randomForest_tune(datasets = list(train_CD[,c(1,6:12)]),
                                 ntree.range = c(50,100,200,300,400,500,600,800,1000),
                                 return.model = T, parallel.cores = -1, seed = 1)

mod_nc_dist <- randomForest_tune(datasets = list(train_NC[,c(1,6:12)]),
                                 ntree.range = c(50,100,200,300,400,500,600,800,1000),
                                 return.model = T, parallel.cores = -1, seed = 1)

save(mod_cd_dist, mod_nc_dist, file = "mod_cd_nc_dist.RData")

library(randomForest)
library(pROC)

res_cd_dist <- predict(mod_cd_dist, test_CD, type = "prob")
res_nc_dist <- predict(mod_nc_dist, test_NC, type = "prob")

roc_cd_dist <- pROC::roc(test_CD$label, res_cd_dist[,1])
roc_nc_dist <- pROC::roc(test_NC$label, res_nc_dist[,1])


for (i in ls(pattern = "res_cd")) {
        res <- get(i)
        res_data <- pROC::roc(test_CD$label,  res[,1])
        assign(paste0("roc_", gsub("res_", "", i)), res_data)
        print(paste0(i, ": ", round(pROC::auc(test_CD$label,  res[,1]), 4)))
}

for (i in ls(pattern = "res_nc")) {
        res <- get(i)
        res_data <- pROC::roc(test_NC$label,  res[,1])
        assign(paste0("roc_", gsub("res_", "", i)), res_data)
        print(paste0(i, ": ", round(pROC::auc(test_NC$label,  res[,1]), 4)))
}

save(list = c(ls(pattern = "roc_"), ls(pattern = "res_"), "test_CD", "test_NC", "train_CD", "train_NC"),
     file = "ROC_res_v1.7.1.RData")

mod_nc_consv_dist <- randomForest_tune(datasets = list(train_NC[,c(1:12)]),
                                      ntree.range = c(50,100,200,300,400,500,600,800,1000),
                                      return.model = T, parallel.cores = -1, seed = 1)

mod_nc_consv_dist_vep <- randomForest_tune(datasets = list(train_NC),
                                          ntree.range = c(50,100,200,300,400,500,600,800,1000),
                                          return.model = T, parallel.cores = -1, seed = 1)


mod_cd_consv_dist_vepAll <- randomForest_tune(datasets = list(train_CD),
                                      ntree.range = c(50,100,200,300,400,500,600,800,1000),
                                      return.model = T, parallel.cores = -1, seed = 1)

mod_cd_consv_dist_vepCons <- randomForest_tune(datasets = list(train_CD[,c(1:49)]),
                                             ntree.range = c(50,100,200,300,400,500,600,800,1000),
                                             return.model = T, parallel.cores = -1, seed = 1)

mod_cd_consv_dist <- randomForest_tune(datasets = list(train_CD[,c(1:12)]),
                                      ntree.range = c(50,100,200,300,400,500,600,800,1000),
                                      return.model = T, parallel.cores = -1, seed = 1)

save(list = ls(pattern = "mod_"), file = "mod_cd_nc_dist.RData")

res_cd_consv_dist_vepAll <- predict(mod_cd_consv_dist_vepAll, test_CD, type = "prob")
res_cd_consv_dist_vepCons <- predict(mod_cd_consv_dist_vepCons, test_CD, type = "prob")
res_cd_consv_dist <- predict(mod_cd_consv_dist, test_CD, type = "prob")

res_nc_consv_dist <- predict(mod_nc_consv_dist, test_NC, type = "prob")
res_nc_consv_dist_vep <- predict(mod_nc_consv_dist_vep, test_NC, type = "prob")


pROC::ggroc(data = list(consv = roc_cd_consv, 
                        consv_vepAll = roc_cd_consv_vepAll,
                        consv_vepCons = roc_cd_consv_vepCons,
                        consv_dist = roc_cd_consv_dist,
                        consv_dist_vepAll = roc_cd_consv_dist_vepAll,
                        consv_dist_vepCons = roc_cd_consv_dist_vepCons,
                        vepAA = roc_cd_vep_aa, 
                        vepAll = roc_cd_vep_all, 
                        vepCons = roc_cd_vep_cons,
                        dist = roc_cd_dist), 
            legacy.axes = T, size = 1)

pROC::ggroc(data = list(consv = roc_nc_consv,
                        consv_vepCons = roc_nc_consv_vepCons,
                        consv_dist = roc_nc_consv_dist,
                        consv_dist_vepCons = roc_nc_consv_dist_vep,
                        vepCons = roc_nc_vep_cons,
                        dist = roc_nc_dist), 
            legacy.axes = T, size = 1)  # 6 * 8

write.csv(mod_cd_consv_dist_vepCons$importance, file = "importanceScore.csv")

###### Spectrum ######

chrX_seq <- seqinr::read.fasta("chrX.fa", seqonly = T)
chrX_seq <- unlist(chrX_seq)
chrX_seq <- seqinr::s2c(chrX_seq)
chrX_seq <- toupper(chrX_seq)

check_POS <- function(pos, ref_allele, genome) {
        ref_genome <- genome[pos]
        if (toupper(ref_genome) == toupper(ref_allele)) {
                res <- T
        } else {
                res <- F
                message("Input: ", ref_allele, " Genome: ", ref_genome)
        }
        res
}

load("~/UoB_MSc_Proj/dataset_v1.7_consv_dist_vep_sampled.RData")
vcf <- data_nc

# res_chrX_pos <- function(vcf, genome) {
#         num_variant <- nrow(vcf)
#         res <- sapply(1:num_variant, function(idx, chrX_seq, num_variant) {
#                 message(idx, "/", num_variant)
#                 res <- check_POS(pos = vcf$POS[idx], ref_allele = vcf$REF[idx], genome = chrX_seq)
#         }, chrX_seq = chrX_seq, num_variant = num_variant)
#         # table(res)
# }

compute_kmer <- function(pos, mutation, genome, width, max_k, freq, mode = c("raw", "Euc")) {
        ref_substring <- genome[(pos - width):(pos + width)]
        ref_count <- c()
        out <- c()
        for (wordsize in 1:min(max_k, ((2 * width) + 1))) {
                res_count <- seqinr::count(tolower(ref_substring), wordsize = wordsize, 
                                           freq = freq)
                ref_count <- c(ref_count, res_count)
        }
        
        ref_substring[(width + 1)] <- mutation
        mut_count <- c()
        for (wordsize in 1:min(max_k, ((2 * width) + 1))) {
                res_count <- seqinr::count(toupper(ref_substring), wordsize = wordsize, 
                                           freq = freq, alphabet = seqinr::s2c("ACGT"))
                mut_count <- c(mut_count, res_count)
                if (mode == "Euc") {
                        euc_res <- sqrt(sum((mut_count - ref_count) ^ 2))
                        names(euc_res) <- paste0("Euc_k.", wordsize)
                        out <- c(out, euc_res)
                }
        }
        if (mode == "raw") {
                out <- c(ref_count, mut_count)
        } 
        out
}

compute_spectrum <- function(vcf, genome, input_width, label = NULL, freq, mode) {
        num_variant <- nrow(vcf)
        
        # res <- sapply(1:num_variant, function(idx, genome, num_variant) {
        #         message(idx, "/", num_variant)
        #         res <- check_POS(pos = vcf$POS[idx], ref_allele = vcf$REF[idx], genome = genome)
        # }, genome = genome, num_variant = num_variant)
        # 
        # if(!all(res)) stop("Input VCF cannot match reference Genome!")
        
        count_res <- sapply(1:num_variant, function(idx, genome, num_variant, input_width) {
                # message(idx, "/", num_variant)
                res <- compute_kmer(pos = vcf$POS[idx], mutation = vcf$ALT[idx],
                                    genome = genome, width = input_width, max_k = 4, freq = freq, mode = mode)
        }, genome = genome, num_variant = num_variant, input_width = input_width)
        count_res_df <- data.frame(t(count_res))
        if (!is.null(label)) count_res_df <- cbind(label = label, count_res_df)
        count_res_df
}

format_spectrum <- function(vcf, genome, max_width = 4, label = vcf$label, 
                            freq = T, mode = c("raw", "Euc")) {
        res <- list()
        for (input_width in 0:max_width) {
                message("Width: ", input_width)
                tmp_res <- compute_spectrum(vcf = vcf, genome = genome, 
                                            input_width = input_width, label = label,
                                            freq = freq, mode = mode)
                res <- c(res, list(tmp_res))
        }
        
        names(res) <- paste0("width.", (2 * (0:max_width) + 1))
        res
}

cd_spect_freq_euc <- format_spectrum(vcf = data_cd, genome = chrX_seq, max_width = 4, 
                                 label = data_cd$label, freq = T, mode = "Euc")
nc_spect_freq_euc <- format_spectrum(vcf = data_nc, genome = chrX_seq, max_width = 4, 
                                 label = data_nc$label, freq = T, mode = "Euc")

cd_spect_num_euc <- format_spectrum(vcf = data_cd, genome = chrX_seq, max_width = 4, 
                                     label = data_cd$label, freq = F, mode = "Euc")
nc_spect_num_euc <- format_spectrum(vcf = data_nc, genome = chrX_seq, max_width = 4, 
                                     label = data_nc$label, freq = F, mode = "Euc")

save(cd_spect, nc_spect, file = "spectNum_cd_nc.RData")
save(cd_spect_freq, nc_spect_freq, file = "spectFreq_cd_nc.RData")

compute_oneGC <- function(pos, mutation, genome, width) {
        ref_substring <- genome[(pos - width):(pos + width)]
        ref_gc <- seqinr::GC(ref_substring, NA.GC = 0)
        
        ref_substring[(width + 1)] <- mutation
        mut_gc <- seqinr::GC(ref_substring, NA.GC = 0)
        
        out <- c(ref_gc = ref_gc, alt_gc = mut_gc)
        out
}

compute_gc <- function(vcf, genome, input_width, label = NULL) {
        num_variant <- nrow(vcf)
        
        # res <- sapply(1:num_variant, function(idx, genome, num_variant) {
        #         message(idx, "/", num_variant)
        #         res <- check_POS(pos = vcf$POS[idx], ref_allele = vcf$REF[idx], genome = genome)
        # }, genome = genome, num_variant = num_variant)
        # 
        # if(!all(res)) stop("Input VCF cannot match reference Genome!")
        
        gc_res <- sapply(1:num_variant, function(idx, genome, num_variant, input_width) {
                # message(idx, "/", num_variant)
                res <- compute_oneGC(pos = vcf$POS[idx], mutation = vcf$ALT[idx],
                                     genome = genome, width = input_width)
        }, genome = genome, num_variant = num_variant, input_width = input_width)
        gc_res_df <- data.frame(t(gc_res))
        if(!is.null(label)) gc_res_df <- cbind(label = label, gc_res_df)
        gc_res_df
}

format_gc <- function(vcf, genome, max_width = 4, label = vcf$label) {
        res <- list()
        for (input_width in 0:max_width) {
                message("Width: ", input_width)
                tmp_res <- compute_gc(vcf = vcf, genome = genome, 
                                      input_width = input_width, label = label)
                res <- c(res, list(tmp_res))
        }
        
        names(res) <- paste0("width.", (2 * (0:max_width) + 1))
        res
}

cd_gc <- format_gc(vcf = data_cd, genome = chrX_seq, max_width = 4, label = data_cd$label)
nc_gc <- format_gc(vcf = data_nc, genome = chrX_seq, max_width = 4, label = data_nc$label)

save(cd_gc, nc_gc, file = "GC_cd_nc.RData")

eval_spectrum <- function(input_spectrum) {
        # browser()
        perf_cv_width <- c()
        for (dataset_idx in 1:length(input_spectrum)) {
                # message(names(input_spectrum)[dataset_idx])
                input_dataset_raw <- input_spectrum[[dataset_idx]]
                k_type <- unique(nchar(names(input_dataset_raw)[-1]))
                perf_cv_k <- list()
                for (dataset_part in 1:length(k_type)) {
                        message(names(input_spectrum)[dataset_idx], "  k.", k_type[1:dataset_part])
                        k_input <- k_type[1:dataset_part]
                        input_dataset <- input_dataset_raw[,c(1, which(nchar(names(input_dataset_raw)) %in% k_input))]
                        perf_cv <- validate_rf(dataset = input_dataset,
                                               ntree.range = c(50,100,200,300,400,500,600,800),
                                               parallel.cores = -1, seed = 1)
                        perf_cv_k <- c(perf_cv_k, list(perf_cv))
                }
                names(perf_cv_k) <- paste0("perf_cv_k.", k_type)
                perf_cv_width <- c(perf_cv_width, list(perf_cv_k))
        }
        names(perf_cv_width) <- paste0("perf_cv_", names(input_spectrum))
        perf_cv_width
}







eval_spectrum_Euc <- function(input_spectrum) {
        # browser()
        perf_cv_width <- c()
        width_type <- c(1,3,5,7,9)
        for (dataset_idx in 1:length(input_spectrum)) {
                # message(names(input_spectrum)[dataset_idx])
                input_dataset_raw <- input_spectrum[[dataset_idx]]
                perf_cv_k <- list()
                k_type <- c(1,2,3,4,5,"all")
                for (dataset_part in c(2,3,4,5,6)) {
                        if (dataset_part == 6) dataset_part <- c(2:5)
                        message(" k.", dataset_part)

                        input_dataset <- input_dataset_raw[,c(1,dataset_part)]
                        perf_cv <- randomForest_tune(datasets = list(input_dataset),
                                                     ntree.range = c(50,100,200,300,400,500),
                                                     return.model = F, parallel.cores = -1, seed = 1)
                        perf_cv_k <- c(perf_cv_k, list(perf_cv))
                }
                names(perf_cv_k) <- paste0("perf_cv_k.", k_type)
                perf_cv_width <- c(perf_cv_width, list(perf_cv_k))
        }
        names(perf_cv_width) <- paste0("perf_cv_width.", width_type)
        perf_cv_width
}

perf_spectrumNum_cd <- eval_spectrum(input_spectrum = cd_spect)
perf_spectrumNum_nc <- eval_spectrum(input_spectrum = nc_spect)

perf_spectrumFreq_cd <- eval_spectrum(input_spectrum = cd_spect_freq)
perf_spectrumFreq_nc <- eval_spectrum(input_spectrum = nc_spect_freq)

perf_spectrumNum_Euc_cd <- eval_spectrum_Euc(input_spectrum = cd_spect_num_euc)
perf_spectrumNum_Euc_nc <- eval_spectrum_Euc(input_spectrum = nc_spect_num_euc)

perf_spectrumFreq_Euc_cd <- eval_spectrum_Euc(input_spectrum = cd_spect_freq_euc)
perf_spectrumFreq_Euc_nc <- eval_spectrum_Euc(input_spectrum = nc_spect_freq_euc)

save(perf_spectrumNum_cd, perf_spectrumNum_nc, perf_spectrumFreq_cd, perf_spectrumFreq_nc,
     file = "perf_spectrum_Num_Freq_nc_cd.RData")

out_perf <- function(input_spectrum_perf, path) {
        perf_spectrumNum_nc <- input_spectrum_perf
        outPerf <- data.frame()
        for(i in 1:length(perf_spectrumNum_nc)) {
                names1 <- gsub("perf_cv_", "", names(perf_spectrumNum_nc[i]))
                for(j in 1:length(perf_spectrumNum_nc[[i]])) {
                        names2 <- gsub("perf_cv_", "", names(perf_spectrumNum_nc[[i]][j]))
                        idx <- row.names(perf_spectrumNum_nc[[i]][[j]]$performance) == paste0("ntree_", perf_spectrumNum_nc[[i]][[j]]$ntree)
                        res <- perf_spectrumNum_nc[[i]][[j]]$performance[idx,,drop = F]
                        res$type <- paste(names1, names2, sep = "_")
                        outPerf <- rbind(outPerf, res)
                }
        }
        outPerf$width <- c(1, 3, 3, 3, 5, 5, 5, 5, 7, 7, 7, 7, 9, 9, 9, 9)
        outPerf$k <- c(1, 1:3, 1:4, 1:4, 1:4)
        # outPerf <- rbind(outPerf, c(0,0,0,0,0,0,1,2), c(0,0,0,0,0,0,1,3), 
        #                  c(0,0,0,0,0,0,1,4), c(0,0,0,0,0,0,3,4))
        outPerf$width <- as.character(outPerf$width)
        outPerf$k <- as.character(outPerf$k)
        # write.csv(outPerf, file = path, quote = F)
        outPerf
}

overall_spectrumFreq_cd <- out_perf(perf_spectrumFreq_cd, path = "overall_spectrumFreq_cd.csv")
overall_spectrumFreq_nc <- out_perf(perf_spectrumFreq_nc, path = "overall_spectrumFreq_nc.csv")
overall_spectrumNum_cd  <- out_perf(perf_spectrumNum_cd,  path = "overall_spectrumNum_cd.csv")
overall_spectrumNum_nc  <- out_perf(perf_spectrumNum_nc,  path = "overall_spectrumNum_nc.csv")

eval_gc <- function(input_gc_dataset) {
        # browser()
        perf_cv_width <- c()
        for (dataset_idx in 1:length(input_gc_dataset)) {
                # message(names(input_spectrum)[dataset_idx])
                input_dataset <- input_gc_dataset[[dataset_idx]]
                width_type <- names(input_gc_dataset[dataset_idx])
                perf_cv <- validate_rf(dataset = input_dataset,
                                       ntree.range = c(50,100,200,300,400,500,600,800),
                                       parallel.cores = -1, seed = 1)
                names(perf_cv) <- paste0("gc_", width_type)
                perf_cv_width <- c(perf_cv_width, list(perf_cv))
        }
        perf_cv_width
}

perf_gc_cd <- eval_gc(input_gc_dataset = cd_gc)
perf_gc_nc <- eval_gc(input_gc_dataset = nc_gc)

out_perf <- function(input_gc_perf, path) {
        perf_spectrumNum_nc <- input_gc_perf
        outPerf <- data.frame()
        for(i in 1:length(perf_spectrumNum_nc)) {
                names1 <- names(perf_spectrumNum_nc[[i]][1])
                idx <- row.names(perf_spectrumNum_nc[[i]][[2]]) == paste0("ntree_", perf_spectrumNum_nc[[i]][[1]])
                res <- perf_spectrumNum_nc[[i]][[2]][idx,,drop = F]
                res$type <- names1
                outPerf <- rbind(outPerf, res)
        }
        write.csv(outPerf, file = path, quote = F)
        outPerf
}

overall_gc_cd <- out_perf(perf_gc_cd, path = "overall_gc_cd.csv")
overall_gc_nc <- out_perf(perf_gc_nc, path = "overall_gc_nc.csv")

library(ggplot2)

freq_cd <- ggplot(data = overall_spectrumFreq_cd, aes(x = width, y = Accuracy, fill = k)) +
        scale_fill_brewer(palette = "Pastel1") +
        geom_bar(stat = "identity", position = position_dodge2(preserve = "single"), color = "grey40") +
        facet_grid(~width, scales = "free_x", space = "free", switch = "x") +
        coord_cartesian(ylim = c(0.4, 0.7)) + 
        scale_y_continuous(breaks = seq(0.4, 0.7, 0.05)) +
        theme(axis.text.x  = element_text(size = 12),
              axis.title.x = element_text(size = 12),
              axis.title.y = element_text(size = 12),
              axis.text.y  = element_text(size = 12),
              legend.title = element_text(size = 12),
              legend.text  = element_text(size = 12),
              legend.position = "bottom")

freq_nc <- ggplot(data = overall_spectrumFreq_nc, aes(x = width, y = Accuracy, fill = k)) +
        scale_fill_brewer(palette = "Pastel1") +
        geom_bar(stat = "identity", position = position_dodge2(preserve = "single"), color = "grey40") +
        facet_grid(~width, scales = "free_x", space = "free", switch = "x") +
        coord_cartesian(ylim = c(0.5, 0.7)) + 
        scale_y_continuous(breaks = seq(0.5, 0.7, 0.05)) +
        theme(axis.text.x  = element_text(size = 12),
              axis.title.x = element_text(size = 12),
              axis.title.y = element_text(size = 12),
              axis.text.y  = element_text(size = 12),
              legend.title = element_text(size = 12),
              legend.text  = element_text(size = 12),
              legend.position = "bottom")

multiplot <- function(..., plotlist = NULL, file, cols = 1, layout = NULL) {
        library(grid)
        
        # Make a list from the ... arguments and plotlist
        plots <- c(list(...), plotlist)
        
        numPlots = length(plots)
        
        # If layout is NULL, then use 'cols' to determine layout
        if (is.null(layout)) {
                # Make the panel
                # ncol: Number of columns of plots
                # nrow: Number of rows needed, calculated from # of cols
                layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                                 ncol = cols, nrow = ceiling(numPlots/cols))
        }
        
        if (numPlots==1) {
                print(plots[[1]])
                
        } else {
                # Set up the page
                grid.newpage()
                pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
                
                # Make each plot, in the correct location
                for (i in 1:numPlots) {
                        # Get the i,j matrix positions of the regions that contain this subplot
                        matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
                        
                        print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                                        layout.pos.col = matchidx$col))
                }
        }
}  

multiplot(freq_cd, freq_nc, cols = 1) # 8 * 10 in

for (i in ls(pattern = "perf_cv_")) {
        write.csv(get(i), file = paste0(i, ".csv"))
}

data_gc <- read.csv("ggplot2_gc_res.csv")

data_gc$Width <- as.factor(data_gc$Width)

ggplot(data = data_gc, aes(x = Width, y = Performance, group = Group, 
                                 colour = Region)) +
        geom_line(size = 1, aes(linetype = Metric)) + 
        coord_cartesian(ylim = c(0.3, 0.65), xlim = c(1.1,4.9)) + 
        scale_y_continuous(breaks = seq(0.3, 0.65, 0.05)) +
        theme(axis.text.x  = element_text(size = 12),
              axis.title.x = element_text(size = 12),
              axis.title.y = element_text(size = 12),
              axis.text.y  = element_text(size = 12),
              legend.title = element_text(size = 12),
              legend.text  = element_text(size = 12),
              legend.position = "bottom") +
        guides(fill = guide_legend(ncol = 2))

###### Distance v2 ######

gencode_chrx <- read.table("C:/Users/alexh/Desktop/TB2/multiz_hg38/gencode.v31.annotation.chrX.gtf", sep = "\t", stringsAsFactors = F)

cal_one_dist <- function(input_type, input_POS, input_df) {
        new_df <- input_df[input_df$V3 == input_type,4:5]
        
        s <- new_df$V4
        e <- new_df$V5
        
        if (any(input_POS <= e & input_POS >= s)) {
                distance <- 1
        } else {
                distance <- 1 / (min(c(abs(input_POS - s), abs(input_POS - e))) + 1)
        }
        distance
}

dist_region <- function(input_POS, gtf_data, threshold) {
        distance_all <- c(CDS = 0, exon = 0, gene = 0, start_codon = 0, stop_codon = 0, transcript = 0, UTR = 0)
        
        s <- gtf_data$V4
        e <- gtf_data$V5
        
        case1 <- which(s >= (input_POS - threshold) & s <= (input_POS + threshold))
        case2 <- which(e >= (input_POS - threshold) & e <= (input_POS + threshold))
        case3 <- which(s <= (input_POS - threshold) & e >= (input_POS + threshold))
        case_all <- unique(c(case1, case2, case3))
        
        input_df <- gtf_data[case_all,]
        
        if (nrow(input_df) != 0) {
                res <- sapply(unique(input_df$V3), cal_one_dist, input_POS = input_POS, 
                              input_df = input_df)
                for (i in 1:length(res)) {
                        distance_all[names(distance_all) == names(res)[i]] <- res[i]
                }
        }
        distance_all
}

cal_distance <- function(all_input_POS, gtf_data = gencode_chrx, threshold, 
                      label = NULL, parallel.cores = 2, cl = NULL) {
        if (is.null(cl)) {
                cl <- parallel::makeCluster(parallel.cores)
                parallel::clusterExport(cl, "cal_one_dist")
        }
        
        distance_test <- parallel::parSapply(cl, all_input_POS, dist_region, 
                                             gtf_data = gtf_data, threshold = threshold)
        
        if (is.null(cl)) {
                parallel::stopCluster(cl)
        }
        
        distance_test_df <- data.frame(t(distance_test), stringsAsFactors = F)
        if (!is.null(label)) distance_test_df <- cbind(label = label, distance_test_df, stringsAsFactors = F)
        distance_test_df
}

save.image("distance_data_helix.RData")

load("distance_data_helix.RData")

rm(list = ls(pattern = "distanceFeatureSet_"))

cl <- parallel::makeCluster(20)
parallel::clusterExport(cl, "cal_one_dist")

for (threshold_input in c(1, 10, 50, 100, 200, 300, 500, 700, 1000, 1500, 2000, 2500, 
                          3000, 3500, seq(4000, 10000, 1000))) {
        message("Proccessing: ", threshold_input, "     ", Sys.time())
        
        distanceDataset <- cal_distance(all_input_POS = c(data_nc$POS), 
                                        gtf_data = gencode_chrx, threshold = threshold_input, cl = cl, 
                                        label = c(as.character(data_nc$label)))
        distanceDataset$label <- as.factor(distanceDataset$label)
        assign(paste0("distanceFeatureSet_", threshold_input), distanceDataset)
        save(list = ls(pattern = "distanceFeatureSet_"), file = "distanceFeatureSet_nc.RData")
}
parallel::stopCluster(cl)

#

rm(list = ls(pattern = "distanceFeatureSet_"))

evaluateDistance <- function(mode = c("cd", "nc")) {
        if (mode == "nc") {
                load("distanceFeatureSet_nc.RData")
                
                for (data_name in ls(pattern = "distanceFeatureSet_")) {
                        message("Proccessing: ", data_name, "     ", Sys.time())
                        dataset_input <- get(data_name)
                        res_distance <- validate_rf(dataset = dataset_input, label.col = 1,
                                                    positive.class = "Positive", folds.num = 10,
                                                    ntree.range = c(50, 100, 200, 300, 500, 700, 1000, 1500),
                                                    seed = 1, parallel.cores = 10)
                        assign(paste0("res_", data_name), res_distance)
                        save(list = ls(pattern = "res_distanceFeatureSet_"), file = "res_distance_nc.RData")
                }
                
                perf_distance <- c()
                for (i in ls(pattern = "res_distance", envir = environment())) {
                        res_data <- get(i)
                        perf_out <- res_data$perf
                        perf_out$WindowSize <- as.numeric(gsub("res_distanceFeatureSet_", "", i))
                        perf_distance <- rbind(perf_distance, perf_out)
                }
                
                perf_distance$Region = "NonCoding"
                
                row.names(perf_distance) <- NULL
                perf_distance_nc <- perf_distance
                save(perf_distance_nc, file = "perf_distance_nc.RData")
        } else {
                load("distanceFeatureSet_cd.RData")
                
                for (data_name in ls(pattern = "distanceFeatureSet_")) {
                        message("Proccessing: ", data_name, "     ", Sys.time())
                        dataset_input <- get(data_name)
                        res_distance <- validate_rf(dataset = dataset_input, label.col = 1,
                                                    positive.class = "Positive", folds.num = 10,
                                                    ntree.range = c(50, 100, 200, 300, 500, 700, 1000, 1500),
                                                    seed = 1, parallel.cores = 10)
                        assign(paste0("res_", data_name), res_distance)
                        save(list = ls(pattern = "res_distanceFeatureSet_"), file = "res_distance_cd.RData")
                }
                
                perf_distance <- c()
                for (i in ls(pattern = "res_distance", envir = environment())) {
                        res_data <- get(i)
                        perf_out <- res_data$perf
                        perf_out$WindowSize <- as.numeric(gsub("res_distanceFeatureSet_", "", i))
                        perf_distance <- rbind(perf_distance, perf_out)
                }
                
                perf_distance$Region = "Coding"
                
                row.names(perf_distance) <- NULL
                perf_distance_cd <- perf_distance
                save(perf_distance_cd, file = "perf_distance_cd.RData")
        }
        perf_distance
}

perf_dist_cd <- evaluateDistance("cd")
perf_dist_nc <- evaluateDistance("nc")


test_cd_vcf <- paste(test_cd$CHROM, test_cd$POS, ".", test_cd$REF, test_cd$ALT, sep = "\t")