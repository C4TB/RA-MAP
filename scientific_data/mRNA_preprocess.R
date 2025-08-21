# Set working directory
setwd('~/Documents/QMUL/RA-MAP')

# Load library
library(limma)

# Import data
x1 <- read.ilmn(files = './Data/Raw/baseline_sampleprobeprofile.txt', 
                ctrlfiles = './Data/Raw/baseline_controlprobeprofile.txt',
                other.columns = 'Detection')
x2 <- read.ilmn(files = './Data/Raw/6m_sample_probe.txt', 
                ctrlfiles = './Data/Raw/6m_control_probe.txt',
                other.columns = 'Detection')
x <- merge(x1, x2)

# Normalize
y <- neqc(x)

# Filter
keep <- rowSums(y$other$Detection < 0.05) >= (ncol(y) / 10) # =219
y <- y[keep, ]

# Aggregate
y <- avereps(y, ID = y$genes$SYMBOL)

# Export
saveRDS(y$E, './Data/microarray_dat.rds')
