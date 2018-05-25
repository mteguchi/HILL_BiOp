#fix_UWindex.R

# Originally written for DcVsDGN manuscript. Converted from there.

ifelse(Sys.info()[1] == 'Linux',
       source('~/Documents/R/tools/TomosFunctions.R'),
       source('~/R/tools/TomosFunctions.R'))
library(dplyr)

UWidx.file <- paste0('data/UW36N_19670101_20170130_Daily_data.txt')
UWidx.file2 <- paste0('data/UW36N_19670101_20170130_Daily_data_int.csv')

UWidx.data <- read.table(UWidx.file, header = FALSE)
colnames(UWidx.data) <- c("Date1", "UW_idx")

UWidx.data[UWidx.data$UW_idx == -9999, "UW_idx"] <- NA

UW.int <- vector(mode = "numeric", length = dim(UWidx.data)[1])
# or just linear interpolate
i <- 1
while (i <= dim(UWidx.data)[1]){
  if (!is.na(UWidx.data$UW[i])){
    UW.int[i] <- UWidx.data$UW[i]
    i <- i + 1
  } else {
    non.na.val1 <- UWidx.data$UW[i-1]
    i1 <- i
    i <- i + 1
    while (is.na(UWidx.data$UW[i])) {
      i <- i + 1
    }
    non.na.val2 <- UWidx.data$UW[i]
    i2 <- i - 1

    UW.int[i1:i2] <- non.na.val1 + (((non.na.val2 - non.na.val1)/(i2-i1+2)) * c(1:(i2-i1+1)))
  }
}

UWidx.data$UWint <- UW.int
UWidx.data %>% mutate(year = yyyymmdd2Y(Date1),
                      month = yyyymmdd2M(Date1),
                      date = yyyymmdd2D(Date1)) %>%
  group_by(year) %>%
  mutate(cumulativeUW = cumsum(UWint)) -> UWidx.data

write.table(UWidx.data, file = UWidx.file2,
            quote = F, sep = ",", row.names = F)
