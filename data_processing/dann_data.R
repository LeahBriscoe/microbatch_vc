local_dir = "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/"

meta_1 = readRDS(paste0(local_dir,"Hispanic_k7/metadata.rds"))


meta_2 = readRDS(paste0(local_dir,"AGP_max_k7/metadata.rds"))
require(dplyr)
dim(meta_2)
dim(meta_1)
colnames(meta_1)
table(meta_1$placeofbirth_group.x)
table(meta_1$diabetes_lab_v2.y)
table(meta_2$diabetes)
table(meta_2 %>% filter(body_site.x == "UBERON:feces") %>% select(country_of_birth.x))


# smoker, hypertension, cancer,cholesterol, yogurt education