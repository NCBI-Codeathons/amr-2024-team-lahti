# There will be two analyses in this script.

## Analysis 1. Explore total ARG abundance in adult stool metagenomes determined based on assembly

## Analysis 2. Explore total ARG abundance in adult stool metagenomes determined based on raw reads


# Analysis 1. Explore total ARG abundance in adult stool metagenomes based on assembly
## - You will compare total ARG abundance of stool metagenomes across countries using only the healthy subjects currently not using antibiotics
## - Than you will calculate correlation between the median total ARG abundance and the antibiotic usage of the countries 

## read sample metadata table, ARG abundance matrix, and ARG metadata table.
metadata <- read.table(metadata_file, sep="\t", header=TRUE)
arg_assembly_abund <- read.table(arg_assembly_abund_file, sep="\t", header=TRUE, row.names=1, check.names=FALSE)
arg_info <- read.table(arg_info_file, sep="\t", header=TRUE)

## add the [STUDY__SAMPLE] formatted column to the metadata table; so that sample names in the ARG abundance table can be directly matched to this column.
metadata$study_sample <- paste0(metadata$Study, "__", metadata$Sample_assembly_ID)

## from the metadata table, select samples that pass all of the following criteria
### - all 40 SCGs recovered from the metagenome
### - stool from adult subjects,
### - subjects that are healthy and not taking antibiotics currently,
## you can try to select such samples based on inspecting the 
### - column 'Tier_2.Recover_all_40_SCGs'
### - column 'AgeCategory' and 'BodySite'
### - column 'Disease' and 'antibiotic_current_use_binary'
selected_metadata <- metadata[metadata$Tier_2.Recover_all_40_SCGs == "Yes",]
selected_metadata <- selected_metadata[!is.na(selected_metadata$AgeCategory),]
selected_metadata <- selected_metadata[!is.na(selected_metadata$BodySite),]
selected_metadata <- selected_metadata[!is.na(selected_metadata$Disease),]
selected_metadata <- selected_metadata[!is.na(selected_metadata$antibiotic_current_use_binary),]
selected_metadata <- selected_metadata[selected_metadata$AgeCategory %in% c("School age", "Adult", "Senior"),]
selected_metadata <- selected_metadata[selected_metadata$BodySite == "Stool",]
selected_metadata <- selected_metadata[selected_metadata$Disease == "Healthy",]
selected_metadata <- selected_metadata[selected_metadata$antibiotic_current_use_binary == "no",]
nrow(selected_metadata)
count(selected_metadata$Country)
metadata <- selected_metadata


## force the metadata table & ARG abundance table to have the same number of rows (= samples) with the same order of rows.
nrow(metadata)
nrow(arg_assembly_abund)
ncol(arg_assembly_abund)
sum(metadata$study_sample %in% row.names(arg_assembly_abund))
sum(row.names(arg_assembly_abund) %in% metadata$study_sample)
metadata <- metadata[metadata$study_sample %in% row.names(arg_assembly_abund),]
arg_assembly_abund <- arg_assembly_abund[row.names(arg_assembly_abund) %in% metadata$study_sample,]
nrow(metadata)
nrow(arg_assembly_abund)
metadata <- metadata[order(metadata$study_sample),]
arg_assembly_abund <- arg_assembly_abund[order(row.names(arg_assembly_abund)),]
sum(metadata$study_sample == row.names(arg_assembly_abund))
### the last three numbers should be same


## add total ARG abundance of each sample, to the metadata table as a new column
metadata$total_arg_cpg <- rowSums(arg_assembly_abund)

## For each country, calculate median and inter quartile ranges of the abundance
## also exclude the poorly sampled country here (n < 10)
country_data <- as.data.frame(count(metadata$Country))
colnames(country_data) <- c("country", "n")
country_data$low_quar <- tapply(metadata$total_arg_cpg, metadata$Country, quantile, p=0.25)
country_data$median <- tapply(metadata$total_arg_cpg, metadata$Country, "median")
country_data$up_quar <- tapply(metadata$total_arg_cpg, metadata$Country, quantile, p=0.75)
country_data <- country_data[country_data$n >= 10,]

## Order the countries by descending order by median ARG abundance 
country_data <- country_data[order(-country_data$median),]
ranked_country <- as.vector(country_data$country)

## Adjust metadata table: exclude the samples from the countries that were dropped out, Attach country rank index to the samples.
metadata <- metadata[metadata$Country %in% ranked_country,]
metadata$country_rank <- match(metadata$Country, ranked_country)


## Plot the country comparison data using a boxplot
### color the countries according to the continent
### "CHN" "ITA" "MDG" "FRA" "USA" "KAZ" "ESP" "DEU" "ISR" "SWE" "DNK" "CAN" "AUT" "NLD"
country_color <- c("firebrick", "mediumblue", "darkorange", "mediumblue", "darkkhaki", "firebrick", "mediumblue", "mediumblue", "firebrick", "mediumblue", "mediumblue", "darkkhaki", "mediumblue","mediumblue")
gp <- ggplot(metadata) + theme_bw() + theme(panel.grid=element_blank(), panel.background = element_blank(), panel.border=element_blank(), axis.line = element_line(color="black"))
gp <- gp + geom_abline(slope=0, intercept=0, color="gray80", linetype=2, size=0.2)
gp <- gp + geom_abline(slope=0, intercept=2, color="gray80", linetype=2, size=0.2)
gp <- gp + geom_abline(slope=0, intercept=4, color="gray80", linetype=2, size=0.2)
gp <- gp + geom_abline(slope=0, intercept=6, color="gray80", linetype=2, size=0.2)
gp <- gp + geom_abline(slope=0, intercept=8, color="gray80", linetype=2, size=0.2)
gp <- gp + geom_jitter(aes(x=as.factor(country_rank), y=total_arg_cpg, color=as.factor(country_rank)), size=0.05, alpha=0.2, shape=16, width=0.16) 
gp <- gp + geom_boxplot(aes(x=as.factor(country_rank), y=total_arg_cpg, fill=as.factor(country_rank)), color="gray20", alpha=0.4, width=0.22, outlier.shape=NA, lwd=0.25) 
gp <- gp + scale_color_manual(values=country_color, name="Country")
gp <- gp + scale_fill_manual(values=country_color, name="Country")
gp <- gp + theme(legend.position="none", axis.title=element_text(size=6.5, color="black"), axis.text.x=element_text(size=6.5, color="black", angle=90, hjust=0.9, vjust=0.5), axis.text.y=element_text(size=6.5, color="black"))
gp <- gp + xlab("\nCountry") + ylab("Normalized abundance (cpg)\nof total ARGs\n")
gp <- gp + scale_x_discrete(labels=ranked_country) 
gp <- gp + scale_y_continuous(breaks=c(0, 2, 4, 6, 8))
gp <- gp + coord_cartesian(ylim=c(0, 8.1))
print(gp)
### Save in tiff file
tiff(filename="ARG_abundance.assembly_based.healthy_no_abx_exposure.country_boxplot.tiff", unit="mm", width=61, height=49, res=1200)
print(gp)
dev.off()

## Country vs country pairiwse comparison: wilcoxon test for difference of medians
n_country <- length(ranked_country)
vec_country_a <- c()
vec_country_b <- c()
vec_med_a <- c()
vec_med_b <- c()
vec_fold <- c()
vec_p <- c()
for(ai in 1:(n_country - 1)){
	for(bi in (ai+1):n_country){
		country_a <- ranked_country[ai]
		country_b <- ranked_country[bi]
		a_d <- metadata[metadata$Country == country_a,]
		b_d <- metadata[metadata$Country == country_b,]
		med_a <- median(a_d$total_arg_cpg)
		med_b <- median(b_d$total_arg_cpg)
		twomeds <- c(med_a, med_b)
		fold_diff <- max(twomeds)/min(twomeds)
		mwtest <- wilcox.test(a_d$total_arg_cpg, b_d$total_arg_cpg)
		vec_country_a <- c(vec_country_a, country_a)
		vec_country_b <- c(vec_country_b, country_b)
		vec_med_a <- c(vec_med_a, med_a)
		vec_med_b <- c(vec_med_b, med_b)
		vec_fold <- c(vec_fold, fold_diff)
		vec_p <- c(vec_p, mwtest$p.value)
	}
}
test_res <- data.frame(country_a = vec_country_a, country_b = vec_country_b, median_a = vec_med_a, median_b = vec_med_b, fold_diff = vec_fold, p = vec_p)
test_res$fdr <- p.adjust(test_res$p, "BH")
### Save in tsv file
write.table(test_res, file="ARG_abundance.assembly_based.healthy_no_abx_exposure.country_pairwise_wilcoxon.tsv", sep="\t", quote=FALSE, row.names=FALSE)


## Test for correlation between country medians vs. country antibiotic consumption rates
### First, read the antibiotic usage data
amu_data <- read.table(country_amu_file, sep="\t", header=TRUE)

### There are two sources of country-wise antibiotic consumption data: CDDEP and WHO 
### Let's use CDDEP version as it covers more countries
cddep_data <- amu_data[!is.na(amu_data$CDDEP.DDD_1000pop_year),]

### Exclude the CDDEP data countries that did not have country level ARG data;
### then attach the country-level ARG abundance median/quartiles to the CDDEP data as new columns
cddep_data <- cddep_data[cddep_data$Country %in% country_data$country,]
cddep_data$ARG_median <- country_data$median[match(cddep_data$Country, country_data$country)]
cddep_data$ARG_low_quar <- country_data$low_quar[match(cddep_data$Country, country_data$country)]
cddep_data$ARG_up_quar <- country_data$up_quar[match(cddep_data$Country, country_data$country)]

### Test for normality 
swt <- shapiro.test(as.vector(cddep_data$ARG_median))
swt$statistic
swt$p.value

### Correlation tests: ARG median ~ Country antibiotic usage 
kend_test <- cor.test(cddep_data$ARG_median, cddep_data$CDDEP.DDD_1000pop_year, method="kendall")
kend_test
pear_test <- cor.test(cddep_data$ARG_median, cddep_data$CDDEP.DDD_1000pop_year, method="pearson")
pear_test

### Plot countries on antibiotic usage vs. median abundance x/y scatter plot
country_color <- c("mediumblue", "darkkhaki", "firebrick", "mediumblue", "mediumblue", "mediumblue", "mediumblue", "mediumblue", "firebrick", "mediumblue", "mediumblue", "darkkhaki")
country_color_b <- c("mediumblue", "darkolivegreen", "firebrick", "mediumblue", "mediumblue", "mediumblue", "mediumblue", "mediumblue", "firebrick", "mediumblue", "mediumblue", "darkolivegreen")
amu_max <- max(cddep_data$CDDEP.DDD_1000pop_year)
gene_max <- max(cddep_data$ARG_up_quar)

gp <- ggplot() + theme_bw() + theme(panel.grid=element_blank(), panel.background = element_blank(), panel.border=element_blank(), axis.line = element_line(color="black"))
gp <- gp + stat_smooth(aes(x=cddep_data$CDDEP.DDD_1000pop_year, y=cddep_data$ARG_median), method="glm", fill="gray50", color="gray50", alpha=0.2)
gp <- gp + geom_linerange(aes(x=cddep_data$CDDEP.DDD_1000pop_year, ymin=cddep_data$ARG_low_quar, ymax=cddep_data$ARG_up_quar, color=cddep_data$Country), alpha=0.3)
gp <- gp + geom_point(aes(x=cddep_data$CDDEP.DDD_1000pop_year, y=cddep_data$ARG_median, fill=cddep_data$Country), size=1, alpha=0.8, shape=23, color="black")
gp <- gp + geom_text_repel(aes(x=cddep_data$CDDEP.DDD_1000pop_year, y=cddep_data$ARG_median, label=cddep_data$Country, color=cddep_data$Country), size=1.6) 
gp <- gp + scale_color_manual(values=country_color_b, name="Country")
gp <- gp + scale_fill_manual(values=country_color, name="Country")
gp <- gp + scale_x_continuous(limits=c(0,16000)) + scale_y_continuous(breaks=c(0, 0.5, 1, 1.5, 2))
gp <- gp + theme(legend.position="none", axis.title=element_text(size=6.5, color="black"), axis.text.x=element_text(size=6.5, color="black"), axis.text.y=element_text(size=6.5, color="black"))
gp <- gp + xlab("\nGross antibiotic consumption\nCDDEP statistic (DDD per 1000)") + ylab("Normalized abundance (cpg)\nof total ARGs\n")
tiff(filename="ARG_abundance.assembly_based.healthy_no_abx_exposure.CDDEP_correlation.tiff", unit="mm", width=61, height=49, res=900)
print(gp)
dev.off()

### From the figure you recognize China as obvious outlier, so repeat the correlation tests
cddep_data_exoutlier <- cddep_data[cddep_data$Country != "CHN",]
kend_test <- cor.test(cddep_data_exoutlier$ARG_median, cddep_data_exoutlier$CDDEP.DDD_1000pop_year, method="kendall")
kend_test
pear_test <- cor.test(cddep_data_exoutlier$ARG_median, cddep_data_exoutlier$CDDEP.DDD_1000pop_year, method="pearson")
pear_test

### Also draw again the correlation scatter plot, just changing the smooth line part (smoothing excluding China)
gp <- ggplot() + theme_bw() + theme(panel.grid=element_blank(), panel.background = element_blank(), panel.border=element_blank(), axis.line = element_line(color="black"))
gp <- gp + stat_smooth(aes(x=cddep_data_exoutlier$CDDEP.DDD_1000pop_year, y=cddep_data_exoutlier$ARG_median), method="glm", fill="gray50", color="gray50", alpha=0.2)
gp <- gp + geom_linerange(aes(x=cddep_data$CDDEP.DDD_1000pop_year, ymin=cddep_data$ARG_low_quar, ymax=cddep_data$ARG_up_quar, color=cddep_data$Country), alpha=0.3)
gp <- gp + geom_point(aes(x=cddep_data$CDDEP.DDD_1000pop_year, y=cddep_data$ARG_median, fill=cddep_data$Country), size=1, alpha=0.8, shape=23, color="black")
gp <- gp + geom_text_repel(aes(x=cddep_data$CDDEP.DDD_1000pop_year, y=cddep_data$ARG_median, label=cddep_data$Country, color=cddep_data$Country), size=1.6) 
gp <- gp + scale_color_manual(values=country_color_b, name="Country")
gp <- gp + scale_fill_manual(values=country_color, name="Country")
gp <- gp + scale_x_continuous(limits=c(0,16000)) + scale_y_continuous(breaks=c(0, 0.5, 1, 1.5, 2))
gp <- gp + theme(legend.position="none", axis.title=element_text(size=6.5, color="black"), axis.text.x=element_text(size=6.5, color="black"), axis.text.y=element_text(size=6.5, color="black"))
gp <- gp + xlab("\nGross antibiotic consumption\nCDDEP statistic (DDD per 1000)") + ylab("Normalized abundance (cpg)\nof total ARGs\n")
tiff(filename="ARG_abundance.assembly_based.healthy_no_abx_exposure.CDDEP_correlation.ex_outlier.tiff", unit="mm", width=61, height=49, res=900)
print(gp)
dev.off()







# Analysis 2. Explore total ARG abundance in adult stool metagenomes based on raw reads
## - You will compare total ARG abundance of stool metagenomes across countries using only the healthy subjects currently not using antibiotics
## - Than you will calculate correlation between the median total ARG abundance and the antibiotic usage of the countries 
## - These are exactly same analyses you did with the assembly-based data

## read sample metadata table, ARG abundance matrix, and ARG metadata table.
metadata <- read.table(metadata_file, sep="\t", header=TRUE)
arg_abund <- read.table(arg_read_abund_file, sep="\t", header=TRUE, row.names=1, check.names=FALSE)
arg_info <- read.table(arg_info_file, sep="\t", header=TRUE)

## from the metadata table, select samples that pass all of the following criteria
### - all 40 SCGs recovered from the metagenome
### - stool from adult subjects,
### - subjects that are healthy and not taking antibiotics currently,
## you can try to select such samples based on inspecting the 
### - column 'Tier_2.Recover_all_40_SCGs'
### - column 'AgeCategory' and 'BodySite'
### - column 'Disease' and 'antibiotic_current_use_binary'
selected_metadata <- metadata[metadata$Tier_2.Recover_all_40_SCGs == "Yes",]
selected_metadata <- selected_metadata[!is.na(selected_metadata$AgeCategory),]
selected_metadata <- selected_metadata[!is.na(selected_metadata$BodySite),]
selected_metadata <- selected_metadata[!is.na(selected_metadata$Disease),]
selected_metadata <- selected_metadata[!is.na(selected_metadata$antibiotic_current_use_binary),]
selected_metadata <- selected_metadata[selected_metadata$AgeCategory %in% c("School age", "Adult", "Senior"),]
selected_metadata <- selected_metadata[selected_metadata$BodySite == "Stool",]
selected_metadata <- selected_metadata[selected_metadata$Disease == "Healthy",]
selected_metadata <- selected_metadata[selected_metadata$antibiotic_current_use_binary == "no",]
nrow(selected_metadata)
count(selected_metadata$Country)
metadata <- selected_metadata


## force the metadata table & ARG abundance table to have the same number of rows (= samples) with the same order of rows.
nrow(metadata)
nrow(arg_abund)
ncol(arg_abund)
sum(metadata$Sample_ID %in% row.names(arg_abund))
sum(row.names(arg_abund) %in% metadata$Sample_ID)
metadata <- metadata[metadata$Sample_ID %in% row.names(arg_abund),]
arg_abund <- arg_abund[row.names(arg_abund) %in% metadata$Sample_ID,]
nrow(metadata)
nrow(arg_abund)
metadata <- metadata[order(metadata$Sample_ID),]
arg_abund <- arg_abund[order(row.names(arg_abund)),]
sum(metadata$Sample_ID == row.names(arg_abund))
### the last three numbers should be same


## add total ARG abundance of each sample, to the metadata table as a new column
metadata$total_arg_cpg <- rowSums(arg_abund)

## For each country, calculate median and inter quartile ranges of the abundance
## also exclude the poorly sampled country here (n < 10)
country_data <- as.data.frame(count(metadata$Country))
colnames(country_data) <- c("country", "n")
country_data$low_quar <- tapply(metadata$total_arg_cpg, metadata$Country, quantile, p=0.25)
country_data$median <- tapply(metadata$total_arg_cpg, metadata$Country, "median")
country_data$up_quar <- tapply(metadata$total_arg_cpg, metadata$Country, quantile, p=0.75)
country_data <- country_data[country_data$n >= 10,]

## Order the countries by descending order by median ARG abundance 
country_data <- country_data[order(-country_data$median),]
ranked_country <- as.vector(country_data$country)

## Adjust metadata table: exclude the samples from the countries that were dropped out, Attach country rank index to the samples.
metadata <- metadata[metadata$Country %in% ranked_country,]
metadata$country_rank <- match(metadata$Country, ranked_country)


## Plot the country comparison data using a boxplot
### color the countries according to the continent
### "ESP" "FRA" "CHN" "KAZ" "ISR" "MDG" "DNK" "USA" "ITA" "DEU" "CAN" "AUT" "SWE" "NLD"
country_color <- c("mediumblue","mediumblue", "firebrick","firebrick","firebrick", "darkorange", "mediumblue", "darkkhaki", "mediumblue","mediumblue", "darkkhaki", "mediumblue","mediumblue","mediumblue")
gp <- ggplot(metadata) + theme_bw() + theme(panel.grid=element_blank(), panel.background = element_blank(), panel.border=element_blank(), axis.line = element_line(color="black"))
gp <- gp + geom_abline(slope=0, intercept=0, color="gray80", linetype=2, size=0.2)
gp <- gp + geom_abline(slope=0, intercept=3, color="gray80", linetype=2, size=0.2)
gp <- gp + geom_abline(slope=0, intercept=6, color="gray80", linetype=2, size=0.2)
gp <- gp + geom_abline(slope=0, intercept=9, color="gray80", linetype=2, size=0.2)
gp <- gp + geom_abline(slope=0, intercept=12, color="gray80", linetype=2, size=0.2)
gp <- gp + geom_jitter(aes(x=as.factor(country_rank), y=total_arg_cpg, color=as.factor(country_rank)), size=0.05, alpha=0.2, shape=16, width=0.16) 
gp <- gp + geom_boxplot(aes(x=as.factor(country_rank), y=total_arg_cpg, fill=as.factor(country_rank)), color="gray20", alpha=0.4, width=0.22, outlier.shape=NA, lwd=0.25) 
gp <- gp + scale_color_manual(values=country_color, name="Country")
gp <- gp + scale_fill_manual(values=country_color, name="Country")
gp <- gp + theme(legend.position="none", axis.title=element_text(size=6.5, color="black"), axis.text.x=element_text(size=6.5, color="black", angle=90, hjust=0.9, vjust=0.5), axis.text.y=element_text(size=6.5, color="black"))
gp <- gp + xlab("\nCountry") + ylab("Normalized abundance (cpg)\nof total ARGs\n")
gp <- gp + scale_x_discrete(labels=ranked_country) 
gp <- gp + scale_y_continuous(breaks=c(0, 3, 6, 9, 12))
gp <- gp + coord_cartesian(ylim=c(0, 15))
print(gp)
### Save in tiff file
tiff(filename="ARG_abundance.read_based.healthy_no_abx_exposure.country_boxplot.tiff", unit="mm", width=61, height=49, res=1200)
print(gp)
dev.off()

## Country vs country pairiwse comparison: wilcoxon test for difference of medians
n_country <- length(ranked_country)
vec_country_a <- c()
vec_country_b <- c()
vec_med_a <- c()
vec_med_b <- c()
vec_fold <- c()
vec_p <- c()
for(ai in 1:(n_country - 1)){
	for(bi in (ai+1):n_country){
		country_a <- ranked_country[ai]
		country_b <- ranked_country[bi]
		a_d <- metadata[metadata$Country == country_a,]
		b_d <- metadata[metadata$Country == country_b,]
		med_a <- median(a_d$total_arg_cpg)
		med_b <- median(b_d$total_arg_cpg)
		twomeds <- c(med_a, med_b)
		fold_diff <- max(twomeds)/min(twomeds)
		mwtest <- wilcox.test(a_d$total_arg_cpg, b_d$total_arg_cpg)
		vec_country_a <- c(vec_country_a, country_a)
		vec_country_b <- c(vec_country_b, country_b)
		vec_med_a <- c(vec_med_a, med_a)
		vec_med_b <- c(vec_med_b, med_b)
		vec_fold <- c(vec_fold, fold_diff)
		vec_p <- c(vec_p, mwtest$p.value)
	}
}
test_res <- data.frame(country_a = vec_country_a, country_b = vec_country_b, median_a = vec_med_a, median_b = vec_med_b, fold_diff = vec_fold, p = vec_p)
test_res$fdr <- p.adjust(test_res$p, "BH")
### Save in tsv file
write.table(test_res, file="ARG_abundance.read_based.healthy_no_abx_exposure.country_pairwise_wilcoxon.tsv", sep="\t", quote=FALSE, row.names=FALSE)


## Test for correlation between country medians vs. country antibiotic consumption rates
### First, read the antibiotic usage data
amu_data <- read.table(country_amu_file, sep="\t", header=TRUE)

### There are two sources of country-wise antibiotic consumption data: CDDEP and WHO 
### Let's use CDDEP version as it covers more countries
cddep_data <- amu_data[!is.na(amu_data$CDDEP.DDD_1000pop_year),]

### Exclude the CDDEP data countries that did not have country level ARG data;
### then attach the country-level ARG abundance median/quartiles to the CDDEP data as new columns
cddep_data <- cddep_data[cddep_data$Country %in% country_data$country,]
cddep_data$ARG_median <- country_data$median[match(cddep_data$Country, country_data$country)]
cddep_data$ARG_low_quar <- country_data$low_quar[match(cddep_data$Country, country_data$country)]
cddep_data$ARG_up_quar <- country_data$up_quar[match(cddep_data$Country, country_data$country)]

### Test for normality 
swt <- shapiro.test(as.vector(cddep_data$ARG_median))
swt$statistic
swt$p.value

### Correlation tests: ARG median ~ Country antibiotic usage 
kend_test <- cor.test(cddep_data$ARG_median, cddep_data$CDDEP.DDD_1000pop_year, method="kendall")
kend_test
pear_test <- cor.test(cddep_data$ARG_median, cddep_data$CDDEP.DDD_1000pop_year, method="pearson")
pear_test

### Plot countries on antibiotic usage vs. median abundance x/y scatter plot
country_color <- c("mediumblue", "darkkhaki", "firebrick", "mediumblue", "mediumblue", "mediumblue", "mediumblue", "mediumblue", "firebrick", "mediumblue", "mediumblue", "darkkhaki")
country_color_b <- c("mediumblue", "darkolivegreen", "firebrick", "mediumblue", "mediumblue", "mediumblue", "mediumblue", "mediumblue", "firebrick", "mediumblue", "mediumblue", "darkolivegreen")
amu_max <- max(cddep_data$CDDEP.DDD_1000pop_year)
gene_max <- max(cddep_data$ARG_up_quar)

gp <- ggplot() + theme_bw() + theme(panel.grid=element_blank(), panel.background = element_blank(), panel.border=element_blank(), axis.line = element_line(color="black"))
gp <- gp + stat_smooth(aes(x=cddep_data$CDDEP.DDD_1000pop_year, y=cddep_data$ARG_median), method="glm", fill="gray50", color="gray50", alpha=0.2)
gp <- gp + geom_linerange(aes(x=cddep_data$CDDEP.DDD_1000pop_year, ymin=cddep_data$ARG_low_quar, ymax=cddep_data$ARG_up_quar, color=cddep_data$Country), alpha=0.3)
gp <- gp + geom_point(aes(x=cddep_data$CDDEP.DDD_1000pop_year, y=cddep_data$ARG_median, fill=cddep_data$Country), size=1, alpha=0.8, shape=23, color="black")
gp <- gp + geom_text_repel(aes(x=cddep_data$CDDEP.DDD_1000pop_year, y=cddep_data$ARG_median, label=cddep_data$Country, color=cddep_data$Country), size=1.6) 
gp <- gp + scale_color_manual(values=country_color_b, name="Country")
gp <- gp + scale_fill_manual(values=country_color, name="Country")
gp <- gp + scale_x_continuous(limits=c(0,16000)) + scale_y_continuous(breaks=0:10)
gp <- gp + theme(legend.position="none", axis.title=element_text(size=6.5, color="black"), axis.text.x=element_text(size=6.5, color="black"), axis.text.y=element_text(size=6.5, color="black"))
gp <- gp + xlab("\nGross antibiotic consumption\nCDDEP statistic (DDD per 1000)") + ylab("Normalized abundance (cpg)\nof total ARGs\n")
tiff(filename="ARG_abundance.read_based.healthy_no_abx_exposure.CDDEP_correlation.tiff", unit="mm", width=61, height=49, res=900)
print(gp)
dev.off()

### From the figure you recognize China as obvious outlier (as you did in the assembly-based analysis) so repeat the correlation tests
cddep_data_exoutlier <- cddep_data[cddep_data$Country != "CHN",]
kend_test <- cor.test(cddep_data_exoutlier$ARG_median, cddep_data_exoutlier$CDDEP.DDD_1000pop_year, method="kendall")
kend_test
pear_test <- cor.test(cddep_data_exoutlier$ARG_median, cddep_data_exoutlier$CDDEP.DDD_1000pop_year, method="pearson")
pear_test

### Also draw again the correlation scatter plot, just changing the smooth line part (smoothing excluding China)
gp <- ggplot() + theme_bw() + theme(panel.grid=element_blank(), panel.background = element_blank(), panel.border=element_blank(), axis.line = element_line(color="black"))
gp <- gp + stat_smooth(aes(x=cddep_data_exoutlier$CDDEP.DDD_1000pop_year, y=cddep_data_exoutlier$ARG_median), method="glm", fill="gray50", color="gray50", alpha=0.2)
gp <- gp + geom_linerange(aes(x=cddep_data$CDDEP.DDD_1000pop_year, ymin=cddep_data$ARG_low_quar, ymax=cddep_data$ARG_up_quar, color=cddep_data$Country), alpha=0.3)
gp <- gp + geom_point(aes(x=cddep_data$CDDEP.DDD_1000pop_year, y=cddep_data$ARG_median, fill=cddep_data$Country), size=1, alpha=0.8, shape=23, color="black")
gp <- gp + geom_text_repel(aes(x=cddep_data$CDDEP.DDD_1000pop_year, y=cddep_data$ARG_median, label=cddep_data$Country, color=cddep_data$Country), size=1.6) 
gp <- gp + scale_color_manual(values=country_color_b, name="Country")
gp <- gp + scale_fill_manual(values=country_color, name="Country")
gp <- gp + scale_x_continuous(limits=c(0,16000)) + scale_y_continuous(breaks=0:10)
gp <- gp + theme(legend.position="none", axis.title=element_text(size=6.5, color="black"), axis.text.x=element_text(size=6.5, color="black"), axis.text.y=element_text(size=6.5, color="black"))
gp <- gp + xlab("\nGross antibiotic consumption\nCDDEP statistic (DDD per 1000)") + ylab("Normalized abundance (cpg)\nof total ARGs\n")
tiff(filename="ARG_abundance.read_based.healthy_no_abx_exposure.CDDEP_correlation.ex_outlier.tiff", unit="mm", width=61, height=49, res=900)
print(gp)
dev.off()

