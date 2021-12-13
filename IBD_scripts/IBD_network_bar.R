## Network visualization of IBD sharing between populations with their geo-coordinates on a map
# (inspired by) Ref: https://chrischizinski.github.io/rstats/igraph-ggplotll/

# author: "Dang Liu 25.Nov.2019"

# Last updated: 13.Dec.2021

# Use libraries
library(tidyverse)
library(ggmap)
library(maps)
library(ggrepel)
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

# Read the IBD sharing input and info file
# Modify it according to your path to these files
IBD_path = "/home/dang_liu/Projects/Kula/code_public/IBD_scripts/all.lPSC.Merged" # the merged IBD output from refinedIBD and their merge-ibd-segments.jar
Info_path = "/home/dang_liu/Projects/Kula/code_public/IBD_scripts/Kula.meta.info" # the most important columns are the individual ID (IID), population label (Label), and QC result (Filter)

IBD_raw <- read_delim(IBD_path,delim="\t",col_names=F)
colnames(IBD_raw) <- c("IND1","HAP1","IND2","HAP2","CHR","BEGIN","END","LOD","LEN")
Info <- read_delim(Info_path,delim="\t",col_names=T) %>% filter(Filter=="PASS")
# 1 to 5 cM
IBD_M_L <- IBD_stats(IBD_raw, Info, min=1, max=5)
IBD_M_N <- IBD_stats(IBD_raw, Info, Measure="N", min=1, max=5)
IBD_D_L <- na.omit(melt(IBD_M_L))
colnames(IBD_D_L) <- c("Pop1", "Pop2", "Avg.sum.IBD.L")
IBD_D_N <- na.omit(melt(IBD_M_N))
colnames(IBD_D_N) <- c("Pop1", "Pop2", "Avg.IBD.N")
IBD_1_5 <- IBD_D_L %>% left_join(IBD_D_N) %>% mutate(Range="1to5")
# 5 to 10 cM
IBD_M_L <- IBD_stats(IBD_raw, Info, min=5, max=10)
IBD_M_N <- IBD_stats(IBD_raw, Info, Measure="N", min=5, max=10)
IBD_D_L <- na.omit(melt(IBD_M_L))
colnames(IBD_D_L) <- c("Pop1", "Pop2", "Avg.sum.IBD.L")
IBD_D_N <- na.omit(melt(IBD_M_N))
colnames(IBD_D_N) <- c("Pop1", "Pop2", "Avg.IBD.N")
IBD_5_10 <- IBD_D_L %>% left_join(IBD_D_N) %>% mutate(Range="5to10")
# over 10 cM
IBD_M_L <- IBD_stats(IBD_raw, Info, min=10)
IBD_M_N <- IBD_stats(IBD_raw, Info, Measure="N", min=10)
IBD_D_L <- na.omit(melt(IBD_M_L))
colnames(IBD_D_L) <- c("Pop1", "Pop2", "Avg.sum.IBD.L")
IBD_D_N <- na.omit(melt(IBD_M_N))
colnames(IBD_D_N) <- c("Pop1", "Pop2", "Avg.IBD.N")
IBD_over10 <- IBD_D_L %>% left_join(IBD_D_N) %>% mutate(Range="over10")
# Combine them
IBD <- rbind(IBD_1_5, IBD_5_10, IBD_over10)

# Filter for copies
# Use 0.5 cutoff, so on average at least half of the pairs share IBD blocks
IBD <- IBD %>% filter(Avg.IBD.N > 0.5)


# Define groups
Collingwood_Bay <- c("Northern", "Wanigela", "Airara")
W_Mas <- c("Fergusson", "Normanby", "Mainland_Eastern_Tip")
N_Mas <- c("Trobriand", "Gawa", "Woodlark", "Laughlan")
S_Mas <- c("Misima", "Western_Calvados", "Eastern_Calvados", "Sudest", "Rossel")

Massim <- c(Collingwood_Bay, W_Mas, N_Mas, S_Mas)


# Get the median of sample geo-coordinates for each pop, and make a new info table at population-level
# Filter with lat and long
# For Oceanians
info2 <- Info %>% filter(Label%in%IBD$Pop1 & Label%in%IBD$Pop2 & Region=="Oceania" & is.na(Longitude)==F) %>% group_by(Label) %>% 
  summarise_at(vars(Latitude,Longitude), funs(median(.))) %>% 
  #left_join(select(Info,-(FID:IID),-(Latitude:Longitude))) %>%
  left_join(select(Info,-(IID),-(Latitude:Longitude))) %>%
  distinct(Label, .keep_all = TRUE)
# Make this to keep all the groups shown in the geom_points regardless of the segments filtering below
info3 <- info2


# Combine info and IBD, and make a new data frame with geo-coordinate links
IBD$from.long <- info2$Longitude[match(IBD$Pop1, info2$Label)]
IBD$from.lat <- info2$Latitude[match(IBD$Pop1, info2$Label)]
IBD$to.long <- info2$Longitude[match(IBD$Pop2, info2$Label)]
IBD$to.lat <- info2$Latitude[match(IBD$Pop2, info2$Label)]


# Remove NA segments, within group sharing, and the 0 IBD length sharing pairs
IBD <- IBD %>% drop_na() %>% filter(Pop1!=Pop2) %>% filter(Avg.sum.IBD.L!=0)

# Adjust for the boundary of maximum
IBD$Avg.sum.IBD.L <- ifelse(IBD$Avg.sum.IBD.L >= 30, 30, IBD$Avg.sum.IBD.L)

# Colored by Area
Area = c("PNG_Highland"="#A6CEE3","PNG_Lowland"="#1F78B4","Bismarck_Arch"="#B2DF8A", 
         "Massim"="#33A02C","Solomon_Arch"="#FB9A99","Australia"="#E31A1C","Non_Oceania"="#FDBF6F")
Region = c("Collingwood_Bay"="#8DD3C7","Western_Massim"="#BEBADA","Northern_Massim"="#FB8072","Southern_Massim"="#80B1D3")
# Get the map
map.world <- map_data(map="world")


# Plot!
# The map
p <- ggplot()
p <- p + theme()
p <- p + geom_map(data=map.world, map=map.world, aes(map_id=region), fill="grey", colour="grey", size=0.15)
p <- p + coord_quickmap(ylim=c(min(info3$Latitude),max(info3$Latitude)), xlim=c(min(info3$Longitude),max(info3$Longitude)))

# The egdes
p <- p + geom_segment(data=IBD, aes(x=from.long, xend=to.long, y=from.lat, yend=to.lat, color=Avg.sum.IBD.L), alpha=0.5, size=1)
p <- p + scale_colour_gradientn(colors=c("#4575B4", "#66CCFF", "#FFEDA0", "#D73027"), breaks=c(5,10,15,20,25,30), labels=c("5","10","15","20","25",">=30"))


# New facet label names for range
range.labs <- c("1 to 5 cM", "5 to 10 cM", "Over 10 cM")
names(range.labs) <- c("1to5", "5to10", "over10")
p <- p + facet_wrap(.~Range, nrow=3, labeller = labeller(Range=range.labs))

# The pops
p <- p + geom_point(data=info3, aes(x=Longitude, y=Latitude, fill=Area), pch=21, size=3, alpha=0.75)
p <- p + scale_fill_manual(values=Area)


# The styles
p <- p + theme(panel.background = element_rect(fill = "#333333"), panel.grid = element_blank())
p <- p + theme(legend.text = element_text(size = 12), legend.title = element_text(size = 14))
p <- p + theme(axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14))
p <- p + theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10))
p <- p + theme(strip.text.x = element_text(size = 12), strip.text.y = element_text(size=12))
p <- p + labs(x="Longitude", y="Latitude", fill="Area")
p


#####Update 31.Aug.2021 For the barplots to visualize the quantification of IBD sharing between Massim groups and all the Oceanian groups######
# Groups
Collingwood_Bay <- c("Northern", "Wanigela", "Airara")
W_Mas <- c("Fergusson", "Normanby", "Mainland_Eastern_Tip")
N_Mas <- c("Trobriand", "Gawa", "Woodlark", "Laughlan")
S_Mas <- c("Misima", "Western_Calvados", "Eastern_Calvados", "Sudest", "Rossel")

Massim <- c(Collingwood_Bay, W_Mas, N_Mas, S_Mas)
Outgroup <- c("Africa", "WestEurasia")
Near_Out <- c("EastAsia", "EastAsia_AN", "Australia")
PNG_H <- c("S_Papuan", "Southern_Highlands", "Enga", "Western_Highlands", "Chimbu", "Eastern_Highlands", "Gulf", "Madang_Highland")
PNG_L <- c("Madang_Lowland", "Morobe", "East_Sepik", "Western", "Central")
Bismarck <- c("Manus_New_Ireland", "West_New_Britain", "East_New_Britain")
Solomon <- c("Bougainville", "Vella_Lavella", "Malaita", "Santa_Cruz", "Bellona_Rennell", "Tikopia")

##boxplot for Massim region vs. other groups
IBD <- rbind(IBD_1_5, IBD_5_10, IBD_over10)
IBD <- IBD %>% filter(Avg.IBD.N > 0.5)
info4 <- info3 %>% rename("Pop1"=Label, "Area1"=Area, "Region1"=Region) %>% select(Pop1, Area1, Region1)
info5 <- info3 %>% rename("Pop2"=Label, "Area2"=Area, "Region2"=Region) %>% select(Pop2, Area2, Region2)
d <- IBD %>% left_join(info4) %>% left_join(info5) %>% drop_na()

d$Region1 <- "Others"
d[d$Pop1%in%Collingwood_Bay,]$Region1 <- "Collingwood_Bay"
d[d$Pop1%in%W_Mas,]$Region1 <- "Western_Massim"
d[d$Pop1%in%N_Mas,]$Region1 <- "Northern_Massim"
d[d$Pop1%in%S_Mas,]$Region1 <- "Southern_Massim"

d$Region2 <- "Others"
d[d$Pop2%in%Collingwood_Bay,]$Region2 <- "Collingwood_Bay"
d[d$Pop2%in%W_Mas,]$Region2 <- "Western_Massim"
d[d$Pop2%in%N_Mas,]$Region2 <- "Northern_Massim"
d[d$Pop2%in%S_Mas,]$Region2 <- "Southern_Massim"


Region_col = c("Collingwood_Bay"="#8DD3C7","Western_Massim"="#BEBADA","Northern_Massim"="#FB8072","Southern_Massim"="#80B1D3","Others"="#FFFFFF")
Area_col = c("PNG_Highland"="#A6CEE3","PNG_Lowland"="#1F78B4","Bismarck_Arch"="#B2DF8A", 
         "Massim"="#33A02C","Solomon_Arch"="#FB9A99","Australia"="#E31A1C","Non_Oceania"="#FDBF6F")


d$Pop2 <- factor(d$Pop2, levels=c(Massim, PNG_H, PNG_L, Bismarck, Solomon), ordered=T)
d$Region1 <- factor(d$Region1, levels=c("Collingwood_Bay", "Western_Massim", "Northern_Massim", "Southern_Massim", "Others"), ordered=T)
d$Region2 <- factor(d$Region2, levels=c("Collingwood_Bay", "Western_Massim", "Northern_Massim", "Southern_Massim", "Others"), ordered=T)
d$Area2 <- factor(d$Area2, levels=c("Massim", "PNG_Highland", "PNG_Lowland", "Bismarck_Arch", "Solomon_Arch"), ordered=T)

range.labs <- c("1 to 5 cM", "5 to 10 cM", "Over 10 cM")
names(range.labs) <- c("1to5", "5to10", "over10")

##All
d %>% filter(Region1!="Others") %>% group_by(Pop2, Region1, Range) %>% 
  summarise_at(vars(Avg.sum.IBD.L, Avg.IBD.N), funs(median(.))) %>%
  left_join(select(d,-(c(Avg.sum.IBD.L, Avg.IBD.N, Pop1, Region1, Range)))) %>%
  distinct(Pop2, Region1, Range, .keep_all = TRUE) %>% 
  ggplot() + 
  geom_bar(aes(x=Pop2, y=Avg.sum.IBD.L, color=Area2, fill=Region2), stat="identity") +
  scale_color_manual(values = Area_col) +
  scale_fill_manual(values = Region_col) +
  facet_wrap(.~Region1*Range, nrow=4, ncol=3, labeller = labeller(Range=range.labs)) +
  theme(axis.line.x = element_line(color="black", size = 0.5, linetype = 1),
        axis.line.y = element_line(color="black", size = 0.5, linetype = 1)) +
  theme_bw() +
  #theme(panel.background = element_blank()) +
  labs(x="Group", y="Avg.sum.IBD.L (cM)", fill="Massim region", color="Region") +
  guides(color=guide_legend(override.aes=list(fill="transparent"))) +
  theme(legend.text = element_text(size = 12), legend.title = element_text(size = 14)) +
  theme(axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14)) + 
  theme(axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust=1), axis.text.y = element_text(size = 12)) +
  theme(strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12))
