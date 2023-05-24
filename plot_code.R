##############update 2023/5/24
#####code for plot
library(ggalluvial)
library(dplyr)
setwd("path/test/")
df <- read.table("clinical.txt",sep = "\t",row.names = 1,header = T)
head(df)

mycol <- rep(c("#223D6C","#D20A13","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#029149","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767"),2)
CLI_lodes <- to_lodes_form(df[,1:ncol(df)],
                           key = "id",
                           axes = 1:ncol(df),
                           id = "Cohort")
dim(CLI_lodes)

ggplot(CLI_lodes,
       aes(x = id, stratum = stratum, alluvium = Cohort,
           fill = stratum, label = stratum)) +
  scale_x_discrete(expand = c(0, 0)) + 
  geom_flow(width = 3/8,aes.flow = "backward") +
  coord_flip() + 
  geom_stratum(alpha = .9,width = 3/8) +
  geom_text(stat = "stratum", size = 4,color="black") +
  scale_fill_manual(values = mycol) +
  xlab("") + ylab("") +
  theme_bw() + 
  theme(panel.grid =element_blank()) + 
  theme(panel.border = element_blank()) + 
  theme(axis.line.x = element_blank(),axis.ticks = element_blank(),axis.text.x = element_blank()) + 
  ggtitle("") +
  guides(fill = FALSE) 

ggsave("/path/sankey_stratum_2.pdf")