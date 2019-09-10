
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

for (threshold_input in unique(c(1, 10, 50, 100, 200, 300, 500, 700, 1000, 1500, 2000, 2500, 
                          3000, 3500, seq(4000, 50000, 1000), seq(55000, 100000, 5000)))) {
        message("Proccessing: ", threshold_input, "     ", Sys.time())
        
        distanceDataset <- cal_distance(all_input_POS = c(data_nc$POS), 
                                        gtf_data = gencode_chrx, threshold = threshold_input, cl = cl, 
                                        label = c(as.character(data_nc$label)))
        distanceDataset$label <- as.factor(distanceDataset$label)
        assign(paste0("distanceFeatureSet_", threshold_input), distanceDataset)
        save(list = ls(pattern = "distanceFeatureSet_"), file = "distanceFeatureSet_nc.RData")
}

rm(list = ls(pattern = "distanceFeatureSet_"))

for (threshold_input in unique(c(1, 10, 50, 100, 200, 300, 500, 700, 1000, 1500, 2000, 2500, 
                                 3000, 3500, seq(4000, 50000, 1000), seq(55000, 100000, 5000)))) {
        message("Proccessing: ", threshold_input, "     ", Sys.time())
        
        distanceDataset <- cal_distance(all_input_POS = c(data_cd$POS), 
                                        gtf_data = gencode_chrx, threshold = threshold_input, cl = cl, 
                                        label = c(as.character(data_cd$label)))
        distanceDataset$label <- as.factor(distanceDataset$label)
        assign(paste0("distanceFeatureSet_", threshold_input), distanceDataset)
        save(list = ls(pattern = "distanceFeatureSet_"), file = "distanceFeatureSet_cd.RData")
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

perf_distance_cd <- evaluateDistance("cd")
perf_distance_nc <- evaluateDistance("nc")

write.csv(perf_distance_cd, file = "perf_distance_cd.csv")
write.csv(perf_distance_nc, file = "perf_distance_nc.csv")

#

library(ggplot2)

perf_distance <- read.table("res_distance_all_ggplot.csv", sep = ",", header = T)

perf_distance$Window_Size <- as.factor(perf_distance$Window_Size)

ggplot(data = perf_distance, aes(x = Window_Size, y = Performance, group = Group, 
                                 colour = Region)) +
        geom_line(size = 1, aes(linetype = Metric)) + 
        scale_y_continuous(breaks = seq(0, 1, 0.1)) +
        coord_cartesian(ylim = c(0.2, 0.95)) + 
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8),
              axis.title.x = element_text(size = 12),
              axis.title.y = element_text(size = 12),
              axis.text.y  = element_text(size = 12),
              legend.title = element_text(size = 12),
              legend.text  = element_text(size = 12),
              legend.position = "bottom") # 3 * 12 in

###

cd_dist <- distanceFeatureSet_85000
nc_dist <- distanceFeatureSet_95000



