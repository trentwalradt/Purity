################################################################################################
### twalradt August 9, 2016
################################################################################################

library(Flow)
pairs = readRDS(file = "/gpfs/commons/groups/imielinski_lab/projects/Chantal/db/pairs.rds")
job = Job('/gpfs/commons/groups/imielinski_lab/git/mskilab/flows/tasks/PurityPloidy.task', pairs)
