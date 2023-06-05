# Load libaries
library("data.table"); packageVersion("data.table")
library("dplyr"); packageVersion("dplyr")
library("ggplot2"); packageVersion("ggplot2")

# Restore meta data from RDS file
df <- readRDS("../../data/meta_clean_05172023.rds")

str(df)

# Format data.frame for analysis
df$CowId_Farm <- paste0(df$CowId, "_", df$FarmId)
df$CowId_Farm <- as.factor(df$CowId_Farm)

df$Collection_Date <- as.Date(df$Collection_Date)
df$Collection_Year <- format(as.Date(df$Collection_Date, format="%d/%m/%Y"),"%Y")
df$Collection_Month <- format(as.Date(df$Collection_Date, format="%d/%m/%Y"),"%m")

df$Calving_Date <- as.Date(df$CalvingDate)
df$Calving_Year <- format(as.Date(df$Calving_Date, format="%d/%m/%Y"),"%Y")
df$Calving_Month <- format(as.Date(df$Calving_Date, format="%d/%m/%Y"),"%m")

# Subset data.frame to only look at true samples (i.e., those collected from animals)
FARMS = c("Farm A", "Farm B", "Farm C", "Farm D", "Farm E")
df <- df %>% filter(Type == "Sample" & FarmId %in% FARMS)

# Summarize number of animals by farm
df %>%
  group_by(FarmId) %>% 
  summarise(No.Animals = uniqueN(CowId_Farm))

# Summarize number of samples collected from each farm
df %>%
  group_by(FarmId) %>% 
  summarise(No.Samples = n())

# Summarize number of samples by farm and sampling period
df %>%
  group_by(FarmId, Period) %>% 
  summarise(No.Samples = n())

# Summarize mean number of samples by farm
df %>%
  group_by(FarmId, CowId_Farm) %>%
  summarise(n = n()) %>%
  summarise(
    Mean = median(n, na.rm = TRUE),
    stdDev = sd(n, na.rm = TRUE)
  )

# Summarise median number of samples collected from each animal
df %>%
  group_by(CowId_Farm) %>%
  summarise(n = n()) %>%
  summarise(
    Median = median(n, na.rm = TRUE),
    stdDev = sd(n, na.rm = TRUE),
    Min = min(n, na.rm = TRUE),
    Max = max(n, na.rm = TRUE)
  )

# Histogram of number of samples collected from each animal
png("eda_samples_per_animal.png")
df %>%
  group_by(CowId_Farm) %>%
  summarise(n = n()) %>%
  ggplot(aes(x = n)) +
  geom_histogram(bins = 16)
dev.off()

df %>% group_by(FarmId, Breed) %>%
  summarise(n = n())

# Histogram of number of samples collected by farm and year
png("eda_collection_year.png")
df %>% group_by(FarmId, Collection_Year) %>%
  summarise(n = n()) %>%
  ggplot(aes(x = Collection_Year, y = n)) +
  geom_bar(stat="identity") +
  facet_wrap(~FarmId)
dev.off()

df %>% group_by(FarmId, Collection_Year) %>%
  summarise(n = n())

# Histogram of number of samples collected by farm and month
png("eda_collection_month.png")
df %>% group_by(FarmId, Collection_Month) %>%
  summarise(n = n()) %>%
  ggplot(aes(x = Collection_Month, y = n)) +
  geom_bar(stat="identity") +
  facet_wrap(~FarmId)
dev.off()

# Histogram of number of samples collected by farm and year
png("eda_calving_year.png")
df %>% group_by(FarmId, Calving_Year) %>%
  summarise(n = n()) %>%
  ggplot(aes(x = Calving_Year, y = n)) +
  geom_bar(stat="identity") +
  facet_wrap(~FarmId)
dev.off()

# Histogram of number of samples collected by farm and month
png("eda_calving_month.png")
df %>% group_by(FarmId, Calving_Month) %>%
  summarise(n = n()) %>%
  ggplot(aes(x = Calving_Month, y = n)) +
  geom_bar(stat="identity") +
  facet_wrap(~FarmId)
dev.off()

# Summarize median number samples collected during prepartum and postpartum
df %>%
  group_by(Period, CowId_Farm) %>%
  summarise(n = n()) %>%
  summarise(
    Median = median(n, na.rm = TRUE),
    stdDev = sd(n, na.rm = TRUE),
    Min = min(n, na.rm = TRUE),
    Max = max(n, na.rm = TRUE)
  )

# Summarize number of samples collected from each animal
df %>%
  group_by(CowId_Farm) %>%
  summarise(n = n()) %>%
  group_by(n) %>%
  count(n)

# Extract information about missing cow ids and calving dates
missingCowId <- df %>% filter(is.na(CowId))
missingCalvingDate <- df %>% filter(is.na(CalvingDate) & !is.na(CowId))

length(unique(missingCalvingDate$CowId_Farm))

# Of the samples with missing calving dates, summarize number of samples collected from each animal
missingCalvingDate %>% 
  group_by(CowId_Farm) %>% 
  summarise(n = n()) %>%
  group_by(n) %>%
  count(n)

# Number of samples with prepartum DIM values less than -56
df %>%
  ungroup() %>%
  filter(DIM < -56) %>%
  summarise(n = n()) %>%
  summarise(sum(n))

# Number of samples with postpartum DIM values greater than 35
df %>%
  ungroup() %>%
  filter(DIM > 35) %>%
  summarise(n = n()) %>%
  summarise(sum(n))