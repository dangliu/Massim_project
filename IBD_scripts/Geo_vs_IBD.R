# Geo vs. IBD

# author: "Dang Liu 05.Sep.2021"

# Last updated: 13.Dec.2021

# Libraries
library(geosphere)
library(tidyverse)
library(reshape2)


# Function
# IBD_stats function
# Read IBD and info files, filtering the IBD result by min, max, and chr
# Can decide whether calculating by L (block length) or N (block number)
IBD_stats <- function(IBD,Info,Measure="L",Chr_exclude="",min=0,max=Inf){
  # Process the IBD table
  IBD <- IBD %>% filter(LEN > min & LEN <= max & CHR!=Chr_exclude)
  # Process the info table 
  IID_Pop_1 <- Info %>% select(IID,Label) %>% rename(IND1="IID",POP1="Label")
  IID_Pop_2 <- Info %>% select(IID,Label) %>% rename(IND2="IID",POP2="Label")
  Pop_Size <- Info %>% select(Label) %>% count(Label)
  # Combine Pop info for ind in IBD table
  IBD_pop <- IBD %>% left_join(IID_Pop_1) %>% left_join(IID_Pop_2)
  
  # Sum up length and copy number by each ind pair
  IBD_pop <- IBD_pop %>% mutate(N=1) %>% group_by(IND1,IND2) %>% mutate(Summed_L = sum(LEN), Summed_N = sum(N)) %>% 
    select(IND1,POP1,IND2,POP2,Summed_L,Summed_N) %>% distinct(.keep_all = T)
  
  # Sum up length and copy number by each pop pair
  IBD_pop <- IBD_pop %>% group_by(POP1,POP2)%>% mutate(Summed_L = sum(Summed_L), Summed_N = sum(Summed_N)) %>% 
    select(POP1,POP2,Summed_L,Summed_N) %>% distinct(.keep_all = T)
  # Get the average L and N for each pop pair considering their sample size pair
  Pop_Size_1 <- Pop_Size %>% rename(POP1="Label",SIZE1="n")
  Pop_Size_2 <- Pop_Size %>% rename(POP2="Label",SIZE2="n")
  
  IBD_pop <- IBD_pop %>% left_join(Pop_Size_1) %>% left_join(Pop_Size_2) %>% mutate(Avg_Sum_L = Summed_L / (SIZE1 * SIZE2), Avg_Sum_N = Summed_N / (SIZE1 * SIZE2))
  # Make IBD matrix L or N for all pop pair beginning with 0
  Pop_Labels <- Info %>% select(Label) %>% distinct(.keep_all = T)
  IBD_M <- matrix(0, nrow=length(Pop_Labels$Label), ncol=length(Pop_Labels$Label))
  rownames(IBD_M) <- Pop_Labels$Label
  colnames(IBD_M) <- Pop_Labels$Label
  
  if (Measure=="L"){
    # Fill in the IBD info
    for (i in (1:nrow(IBD_pop))){
      POP1 <- IBD_pop[i,]$POP1
      POP2 <- IBD_pop[i,]$POP2
      Avg_Sum_L <- IBD_pop[i,]$Avg_Sum_L
      IBD_M[POP1,POP2] = IBD_M[POP1,POP2] + Avg_Sum_L
      IBD_M[POP2,POP1] = IBD_M[POP2,POP1] + Avg_Sum_L
    }
  }
  if (Measure=="N"){
    # Fill in the IBD info
    for (i in (1:nrow(IBD_pop))){
      POP1 <- IBD_pop[i,]$POP1
      POP2 <- IBD_pop[i,]$POP2
      Avg_Sum_N <- IBD_pop[i,]$Avg_Sum_N
      IBD_M[POP1,POP2] = IBD_M[POP1,POP2] + Avg_Sum_N
      IBD_M[POP2,POP1] = IBD_M[POP2,POP1] + Avg_Sum_N
    }
  }
  return(IBD_M)
}
# IBD similarity
S_IBD <- function(IBD_M){
  S_IBD_M <- IBD_M
  for (i in rownames(S_IBD_M)){
    for(j in colnames(S_IBD_M)) {
      S_IBD_M[i,j] <- (2 * IBD_M[i,j]) / (IBD_M[i,i] + IBD_M[j,j]) # similarity
      #S_IBD_M[i,j] <- 1 - (2 * IBD_M[i,j]) / (IBD_M[i,i] + IBD_M[j,j]) # distance
    }
  }
  return(S_IBD_M)
}


## Make geo-distance matrix
# Read info table
# Modify it according to your path to these files
info <- read.table("/home/dang_liu/Projects/Kula/code_public/IBD_scripts/Kula.meta.info", header=T)
head(info)

# Keep good quality Massim sample with geo-coordinates
info2 <- info %>% filter(Filter=="PASS" & is.na(Area)==FALSE & is.na(Longitude)==FALSE) %>% filter(Area=="Massim")

# Use the median of sample coordinates
info2 <- info2 %>% group_by(Label) %>% 
  summarise_at(vars(Latitude,Longitude), funs(median(.))) %>% 
  #left_join(select(info,-(FID:IID),-(Latitude:Longitude))) %>%
  left_join(select(info,-(IID),-(Latitude:Longitude))) %>%
  distinct(Label, .keep_all = TRUE)

## Geo matrix in kilometers
geo_m <- distm(cbind(info2$Longitude, info2$Latitude)*(1/1000))
rownames(geo_m) <- info2$Label
colnames(geo_m) <- info2$Label

# Define groups
Collingwood_Bay <- c("Northern", "Wanigela", "Airara")
W_Mas <- c("Fergusson", "Normanby", "Mainland_Eastern_Tip")
N_Mas <- c("Trobriand", "Gawa", "Woodlark", "Laughlan")
S_Mas <- c("Misima", "Western_Calvados", "Eastern_Calvados", "Sudest", "Rossel")

Massim <- c(Collingwood_Bay, W_Mas, N_Mas, S_Mas)

# 03.Sep.2021.Update Kula groupping: Kula = Kula + Half_Kula
Kula <- c("Mainland_Eastern_Tip","Normanby","Fergusson","Trobriand", "Gawa", "Woodlark")
Half_Kula <- c("Laughlan","Misima")
Kula <- c(Kula, Half_Kula)
Non_Kula <- c("Northern", "Wanigela", "Airara", "Western_Calvados", "Eastern_Calvados", "Sudest", "Rossel")


# Order the matrix
geo_m_order <- geo_m[Massim,Massim]
Massim_geo_m <- geo_m_order[Massim, Massim]
Massim_eff_geo_m <- Massim_geo_m
Kula_geo_m <- geo_m[Kula,Kula]
Half_Kula_geo_m <- geo_m[c(Half_Kula, Kula), c(Half_Kula, Kula)]
for (n in rownames(Kula_geo_m)){
  for (k in colnames(Kula_geo_m)){
    Half_Kula_geo_m[n,k] <- NA
  } 
}

Non_Kula_geo_m <- Massim_geo_m
for (n in rownames(Half_Kula_geo_m)){
  for (k in colnames(Half_Kula_geo_m)){
    Non_Kula_geo_m[n,k] <- NA
  } 
}

#write.table(geo_m_order,"/r1/people/dang_liu/Projects/Kula/IBD/geo_m_order.csv",sep=",")
#write.table(Kula_geo_m,"/r1/people/dang_liu/Projects/Kula/IBD/Kula_geo_m.csv",sep=",")


## Make IBD similarity matrix
# Read the IBD sharing input and info file
# Modify it according to your path to these files
IBD_path = "/home/dang_liu/Projects/Kula/code_public/IBD_scripts/all.lPSC.Merged" # the merged IBD output from refinedIBD and their merge-ibd-segments.jar
Info_path = "/home/dang_liu/Projects/Kula/code_public/IBD_scripts/Kula.meta.info" # the most important columns are the individual ID (IID), population label (Label), and QC result (Filter)
IBD <- read_delim(IBD_path,delim="\t",col_names=F)
colnames(IBD) <- c("IND1","HAP1","IND2","HAP2","CHR","BEGIN","END","LOD","LEN")
Info <- read_delim(Info_path,delim="\t",col_names=T) %>% filter(Filter=="PASS")

# Add range from 1-5 to over 10
# Calculate the normalized similarity value by dividing the between pop sharing with the average of each within pop sharing
# Order the Matrix and excluded those are not used as all populations
IBD_M <- IBD_stats(IBD,Info,min=1,max=5)
IBD_m_1_5 <- S_IBD(IBD_M)[Massim,Massim]
IBD_M <- IBD_stats(IBD,Info,min=2,max=6)
IBD_m_2_6 <- S_IBD(IBD_M)[Massim,Massim]
IBD_M <- IBD_stats(IBD,Info,min=3,max=7)
IBD_m_3_7 <- S_IBD(IBD_M)[Massim,Massim]
IBD_M <- IBD_stats(IBD,Info,min=4,max=8)
IBD_m_4_8 <- S_IBD(IBD_M)[Massim,Massim]
IBD_M <- IBD_stats(IBD,Info,min=5,max=9)
IBD_m_5_9 <- S_IBD(IBD_M)[Massim,Massim]
IBD_M <- IBD_stats(IBD,Info,min=6,max=10)
IBD_m_6_10 <- S_IBD(IBD_M)[Massim,Massim]
IBD_M <- IBD_stats(IBD,Info,min=10)
IBD_m_over10 <- S_IBD(IBD_M)[Massim,Massim]


## Convert both matrix into dataframe and combine them
# don't count the within-pop sharing

# Add range
ibd_d_1_5 <- IBD_m_1_5 %>% melt() %>% na.omit() %>% arrange(., Var1) %>% setNames(c('Pop1', 'Pop2', 'IBD')) %>% mutate(Range="1-5")
ibd_d_2_6 <- IBD_m_2_6 %>% melt() %>% na.omit() %>% arrange(., Var1) %>% setNames(c('Pop1', 'Pop2', 'IBD')) %>% mutate(Range="2-6")
ibd_d_3_7 <- IBD_m_3_7 %>% melt() %>% na.omit() %>% arrange(., Var1) %>% setNames(c('Pop1', 'Pop2', 'IBD')) %>% mutate(Range="3-7")
ibd_d_4_8 <- IBD_m_4_8 %>% melt() %>% na.omit() %>% arrange(., Var1) %>% setNames(c('Pop1', 'Pop2', 'IBD')) %>% mutate(Range="4-8")
ibd_d_5_9 <- IBD_m_5_9 %>% melt() %>% na.omit() %>% arrange(., Var1) %>% setNames(c('Pop1', 'Pop2', 'IBD')) %>% mutate(Range="5-9")
ibd_d_6_10 <- IBD_m_6_10 %>% melt() %>% na.omit() %>% arrange(., Var1) %>% setNames(c('Pop1', 'Pop2', 'IBD')) %>% mutate(Range="6-10")
ibd_d_over10 <- IBD_m_over10 %>% melt() %>% na.omit() %>% arrange(., Var1) %>% setNames(c('Pop1', 'Pop2', 'IBD')) %>% mutate(Range="over10")


geo_d <- Massim_geo_m %>% melt() %>% na.omit() %>% arrange(., Var1) %>% setNames(c('Pop1', 'Pop2', 'Geo')) %>% mutate(Type="Shortest")

# Add IBD range
All_d <- rbind(ibd_d_1_5, ibd_d_2_6, ibd_d_3_7, ibd_d_4_8, ibd_d_5_9, ibd_d_6_10, ibd_d_over10) %>% left_join(rbind(geo_d)) %>%
  mutate(Type2="Non_Kula") %>% 
  mutate(Type2=replace(Type2, Pop1%in%Kula & Pop2%in%Kula, "Kula"))

range.labs <- c("1 to 5 cM", "2 to 6 cM", "3 to 7 cM", "4 to 8 cM", "5 to 9 cM", "6 to 10 cM", "Over 10 cM")
names(range.labs) <- c("1-5", "2-6", "3-7", "4-8", "5-9", "6-10", "over10")

All_d$Type2 <- factor(All_d$Type2, levels=c("Kula", "Non_Kula"), ordered=T)
Group_color = c("Kula"="#6A3884", "Non_Kula"="#C4B286")

## Plot
jitter <- position_jitter(width = 0.025, height = 0.025)
All_d %>% mutate(IBD=ifelse(IBD>1, 1, IBD)) %>% mutate(IBD=ifelse(IBD<0, 0, IBD)) %>% # remove the little noise and make sure IBD similarity is between 0 and 1
  filter(Pop1!=Pop2) %>% # remove within group results
  ggplot(aes(x = Geo, y = IBD, fill=Type2)) +
  geom_point(size=4, pch=21, position=jitter) +
  geom_smooth(aes(fill=Type2), method = lm, show.legend=F, color="black") +
  labs(x = "Geo distance (km)", y = "IBD similarity", fill="Group") +
  scale_fill_manual(values = Group_color) +
  #scale_y_continuous(breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1.0), limits=c(0,1)) +
  facet_wrap(.~Range, nrow=3, labeller = labeller(Range=range.labs)) +
  theme_bw() +
  theme(legend.position = c(.85,.1)) +
  theme(legend.text = element_text(size = 12), legend.title = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14)) +
  theme(strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12)) +
  theme(legend.background = element_rect(color="black", linetype="solid"))


