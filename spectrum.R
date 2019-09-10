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

perf_spectrumNum_cd <- eval_spectrum(input_spectrum = cd_spect)
perf_spectrumNum_nc <- eval_spectrum(input_spectrum = nc_spect)

perf_spectrumFreq_cd <- eval_spectrum(input_spectrum = cd_spect_freq)
perf_spectrumFreq_nc <- eval_spectrum(input_spectrum = nc_spect_freq)

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
                        # idx <- row.names(perf_spectrumNum_nc[[i]][[j]]$performance) == paste0("ntree_", perf_spectrumNum_nc[[i]][[j]]$ntree)
                        res <- perf_spectrumNum_nc[[i]][[j]]$perf
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
        row.names(outPerf) <- NULL
        # write.csv(outPerf, file = path, quote = F, row.names = F)
        outPerf
}

overall_spectrumFreq_cd <- out_perf(perf_spectrumFreq_cd, path = "overall_spectrumFreq_cd.csv")
overall_spectrumFreq_nc <- out_perf(perf_spectrumFreq_nc, path = "overall_spectrumFreq_nc.csv")
overall_spectrumNum_cd  <- out_perf(perf_spectrumNum_cd,  path = "overall_spectrumNum_cd.csv")
overall_spectrumNum_nc  <- out_perf(perf_spectrumNum_nc,  path = "overall_spectrumNum_nc.csv")


library(ggplot2)

freq_cd <- ggplot(data = overall_spectrumFreq_cd, aes(x = width, y = Accuracy, fill = k)) +
        scale_fill_brewer(palette = "Pastel1") +
        geom_bar(stat = "identity", position = position_dodge2(preserve = "single"), color = "grey40") +
        facet_grid(~width, scales = "free_x", space = "free", switch = "x") +
        coord_cartesian(ylim = c(0.45, 0.65)) + 
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
        coord_cartesian(ylim = c(0.55, 0.7)) + 
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

num_cd <- ggplot(data = overall_spectrumNum_cd, aes(x = width, y = Accuracy, fill = k)) +
        scale_fill_brewer(palette = "Pastel1") +
        geom_bar(stat = "identity", position = position_dodge2(preserve = "single"), color = "grey40") +
        facet_grid(~width, scales = "free_x", space = "free", switch = "x") +
        coord_cartesian(ylim = c(0.45, 0.65)) + 
        scale_y_continuous(breaks = seq(0.4, 0.7, 0.05)) +
        theme(axis.text.x  = element_text(size = 12),
              axis.title.x = element_text(size = 12),
              axis.title.y = element_text(size = 12),
              axis.text.y  = element_text(size = 12),
              legend.title = element_text(size = 12),
              legend.text  = element_text(size = 12),
              legend.position = "bottom")

num_nc <- ggplot(data = overall_spectrumNum_nc, aes(x = width, y = Accuracy, fill = k)) +
        scale_fill_brewer(palette = "Pastel1") +
        geom_bar(stat = "identity", position = position_dodge2(preserve = "single"), color = "grey40") +
        facet_grid(~width, scales = "free_x", space = "free", switch = "x") +
        coord_cartesian(ylim = c(0.55, 0.7)) + 
        scale_y_continuous(breaks = seq(0.5, 0.7, 0.05)) +
        theme(axis.text.x  = element_text(size = 12),
              axis.title.x = element_text(size = 12),
              axis.title.y = element_text(size = 12),
              axis.text.y  = element_text(size = 12),
              legend.title = element_text(size = 12),
              legend.text  = element_text(size = 12),
              legend.position = "bottom")

multiplot(num_cd, num_nc, cols = 1) # 8 * 10 in

for (i in ls(pattern = "perf_cv_")) {
        write.csv(get(i), file = paste0(i, ".csv"))
}

###### GC ######

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
        if (!is.null(label)) gc_res_df <- cbind(label = label, gc_res_df)
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
                res <- perf_spectrumNum_nc[[i]][[1]]
                
                res$type <- names1
                outPerf <- rbind(outPerf, res)
        }
        row.names(outPerf) <- NULL
        # write.csv(outPerf, file = path, quote = F, row.names = F)
        outPerf
}

overall_gc_cd <- out_perf(perf_gc_cd, path = "overall_gc_cd.csv")
overall_gc_nc <- out_perf(perf_gc_nc, path = "overall_gc_nc.csv")

data_gc <- read.csv("ggplot2_gc_res.csv")

data_gc$Width <- as.factor(data_gc$Width)

library(ggplot2)

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
        guides(fill = guide_legend(ncol = 2)) # 2 * 10 in

######

cd_spect <- cd_spect$width.3[,c(1:5,86:89)]
nc_spect <- nc_spect$width.9[,c(1:85,342:425)]
