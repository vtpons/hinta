library(dplyr)

ROOT_PATH <- "/path/to/project"

EA4_transmitted <- read.table(paste0(ROOT_PATH, "/output/EA4_transmitted_scores.profile"), header = TRUE)
EA4_nontransmitted <- read.table(paste0(ROOT_PATH, "/output/EA4_nontransmitted_scores.profile"), header = TRUE)
info <- read.table(paste0(ROOT_PATH, "/sample_file.txt"), header = TRUE)

trio_duo <- info[, c("GENOTYPED_ID", "trio_duo")]  # 0 = trio, 1 = duo w/ father, 2 = duo w/ mother
colnames(trio_duo) <- c("IID", "trio_duo")

## TRANSMITTED ##
transmitted <- EA4_transmitted %>%
  mutate(group_id = ceiling(row_number() / 2)) %>%
  group_by(group_id) %>%
  summarise(
    IID = first(IID),
    EA4_PRS_T_FATHER = first(SCORESUM) / 2,
    EA4_PRS_T_MOTHER = last(SCORESUM) / 2,
    EA4_PRS_T_FULL = sum(SCORESUM) / 2
  ) %>%
  select(-group_id)

transmitted$IID <- gsub("1_", "", transmitted$IID)

## NON-TRANSMITTED ##
nontransmitted <- EA4_nontransmitted %>%
  mutate(group_id = ceiling(row_number() / 2)) %>%
  group_by(group_id) %>%
  summarise(
    IID = first(IID),
    EA4_PRS_NT_FATHER = first(SCORESUM) / 2,
    EA4_PRS_NT_MOTHER = last(SCORESUM) / 2
  ) %>%
  select(-group_id)

nontransmitted$IID <- gsub("1_", "", nontransmitted$IID)
nontransmitted <- merge(nontransmitted, trio_duo, by = "IID")
nontransmitted$EA4_PRS_NT_FATHER[nontransmitted$trio_duo == 2] <- NA
nontransmitted$EA4_PRS_NT_MOTHER[nontransmitted$trio_duo == 1] <- NA
nontransmitted$trio_duo <- NULL

## MERGE ##
all_PRS <- merge(transmitted, nontransmitted, by = "IID")
write.table(all_PRS, paste0(ROOT_PATH, "/output/EA4_all_scores.txt"), quote = FALSE, row.names = FALSE)
