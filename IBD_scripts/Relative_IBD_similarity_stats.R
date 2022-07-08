## IBD Hst stats with simulated/shuffled distribution to test

# author: "Dang Liu 30.Aug.2020"

# Last updated: 08.Jul.2022

# Use libraries
library(tidyverse)


### Define functions
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

# Hst (Haplotype statistics, what exactly it does is to calculate relative IBD similarity) function
Hst <- function(IBD_M, all, within, within_exclude=list()){
  # Order the Matrix and excluded those are not used as all populations
  IBD_M <- IBD_M[all,all]
  S_IBD_M <- IBD_M[all,all]
  # Calculate the normalized similarity value by dividing the between pop sharing with the average of each within pop sharing
  for (i in rownames(S_IBD_M)){
    for(j in colnames(S_IBD_M)) {
      S_IBD_M[i,j] <- (2 * IBD_M[i,j]) / (IBD_M[i,i] + IBD_M[j,j])
      #S_IBD_M[i,j] <- (2 * IBD_M[i,j]) / (IBD_M[i,i] + IBD_M[j,j])
    }
  }
  # Hst = (Sw - Sa) / (1 - Sa)
  S_IBD_M_within <- S_IBD_M[within,within]
  if (length(within_exclude)!=0){
    for (r in within_exclude){
      for (c in within_exclude){
        S_IBD_M_within[r,c] <- NA
      }
    }
  }
  S_IBD_M_all <- S_IBD_M
  d_Hst = (mean(S_IBD_M_within[upper.tri(S_IBD_M_within)], na.rm=T) - mean(S_IBD_M_all[upper.tri(S_IBD_M_all)], na.rm=T)) / (1 - mean(S_IBD_M_all[upper.tri(S_IBD_M_all)], na.rm=T))
  return(d_Hst)
}

# Std function
std <- function(x) sd(x)/sqrt(length(x))

# Hst_mean_std function
# By jackknife each of the chr 
Hst_mean_std <- function(IBD, Info, Measure="L", min=0, max=Inf, all, within, within_exclude=list()){
  Hst_list <- sapply(1:22, function(x) Hst(IBD_M=IBD_stats(IBD,Info,Measure=Measure,min=min,max=max,Chr_exclude=x),
                                           all=all,
                                           within=within,
                                           within_exclude=within_exclude))
  Hst_mean = mean(Hst_list)
  Hst_std = std(Hst_list)
  Hst_mean_std <- c(Hst_mean, Hst_std)
  return(Hst_mean_std)
}

# Perform a permutation for randomly assigned the targeted groups from all studied groups, and sample the Hst difference between the target and non-target groups
Sample_Hst_diff <- function(seed, IBD, Info, Measure="L", min=0, max=Inf, all, T_all, within, within_exclude=list()){
  set.seed(seed)
  n = length(within)
  G1 = sample(T_all, n)
  T <- Hst(IBD_M=IBD_stats(IBD,Info,Measure=Measure,min=min,max=max),
           all=all,
           within=G1,
           within_exclude=within_exclude)
  NT <- Hst(IBD_M=IBD_stats(IBD,Info,Measure=Measure,min=min,max=max),
            all=all,
            within=T_all,
            within_exclude=G1)
  diff <- T - NT
  return(diff)
}


# Custom boplot function for determing the quantiles
stat_boxplot_custom <- function(mapping = NULL, data = NULL,
                                geom = "boxplot", position = "dodge",
                                ...,
                                qs = c(.05, .25, 0.5, 0.75, 0.95),
                                na.rm = FALSE,
                                show.legend = NA,
                                inherit.aes = TRUE) {
  layer(
    data = data,
    mapping = mapping,
    stat = StatBoxplotCustom,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      qs = qs,
      ...
    )
  )
}
StatBoxplotCustom <- ggproto("StatBoxplotCustom", Stat,
                             required_aes = c("x", "y"),
                             non_missing_aes = "weight",
                             
                             setup_params = function(data, params) {
                               params$width <- ggplot2:::"%||%"(
                                 params$width, (resolution(data$x) * 0.75)
                               )
                               
                               if (is.double(data$x) && !ggplot2:::has_groups(data) && any(data$x != data$x[1L])) {
                                 warning(
                                   "Continuous x aesthetic -- did you forget aes(group=...)?",
                                   call. = FALSE
                                 )
                               }
                               
                               params
                             },
                             
                             compute_group = function(data, scales, width = NULL, na.rm = FALSE, qs = c(.05, .25, 0.5, 0.75, 0.95)) {
                               
                               if (!is.null(data$weight)) {
                                 mod <- quantreg::rq(y ~ 1, weights = weight, data = data, tau = qs)
                                 stats <- as.numeric(stats::coef(mod))
                               } else {
                                 stats <- as.numeric(stats::quantile(data$y, qs))
                               }
                               names(stats) <- c("ymin", "lower", "middle", "upper", "ymax")
                               iqr <- diff(stats[c(2, 4)])
                               
                               outliers <- (data$y < stats[1]) | (data$y > stats[5])
                               
                               if (length(unique(data$x)) > 1)
                                 width <- diff(range(data$x)) * 0.9
                               
                               df <- as.data.frame(as.list(stats))
                               df$outliers <- list(data$y[outliers])
                               
                               if (is.null(data$weight)) {
                                 n <- sum(!is.na(data$y))
                               } else {
                                 # Sum up weights for non-NA positions of y and weight
                                 n <- sum(data$weight[!is.na(data$y) & !is.na(data$weight)])
                               }
                               
                               df$notchupper <- df$middle + 1.58 * iqr / sqrt(n)
                               df$notchlower <- df$middle - 1.58 * iqr / sqrt(n)
                               
                               df$x <- if (is.factor(data$x)) data$x[1] else mean(range(data$x))
                               df$width <- width
                               df$relvarwidth <- sqrt(n)
                               df
                             }
)




#############

# Read the IBD sharing input and info file
# Modify it according to your path to these files
IBD_path = "/home/dang_liu/Projects/Kula/code_public/IBD_scripts/all.lPSC.Merged" # the merged IBD output from refinedIBD and their merge-ibd-segments.jar
Info_path = "/home/dang_liu/Projects/Kula/code_public/IBD_scripts/Kula.meta.info" # the most important columns are the individual ID (IID), population label (Label), and QC result (Filter)

IBD <- read_delim(IBD_path,delim="\t",col_names=F)
colnames(IBD) <- c("IND1","HAP1","IND2","HAP2","CHR","BEGIN","END","LOD","LEN")
Info <- read_delim(Info_path,delim="\t",col_names=T) %>% filter(Filter=="PASS")

# Define groups
Collingwood_Bay <- c("Northern", "Wanigela", "Airara")
W_Mas <- c("Fergusson", "Normanby", "Mainland_Eastern_Tip")
N_Mas <- c("Trobriand", "Gawa", "Woodlark", "Laughlan")
S_Mas <- c("Misima", "Western_Calvados", "Eastern_Calvados", "Sudest", "Rossel")

Massim <- c(Collingwood_Bay, W_Mas, N_Mas, S_Mas)
Outgroup <- c("Africa", "WestEurasia")
Near_Out <- c("EastAsia", "EastAsia_AN", "Australia")
PNG_H <- c("Southern_Highlands", "Enga", "Western_Highlands", "Chimbu", "Eastern_Highlands", "Madang_Highland")
PNG_L <- c("Madang_Lowland", "Morobe", "East_Sepik", "Western", "Gulf", "Central")
Bismarck <- c("Manus_New_Ireland", "West_New_Britain", "East_New_Britain")
Solomon <- c("Bougainville", "Vella_Lavella", "Malaita", "Santa_Cruz", "Bellona_Rennell", "Tikopia")

Kula <- c("Mainland_Eastern_Tip","Normanby","Fergusson","Trobriand", "Gawa", "Woodlark")
Half_Kula <- c("Laughlan","Misima")
Non_Kula <- c("Northern", "Wanigela", "Airara", "Western_Calvados", "Eastern_Calvados", "Sudest", "Rossel")
All_Kula <- c(Kula, Half_Kula)

# This is the "all" group considering as background
D <- c(PNG_H, PNG_L, Bismarck, Solomon, Massim)

########################

# Assign group to group
assign("Kula", All_Kula)
#assign("Half_Kula", c(Kula, Half_Kula))
assign("Non_Kula", Massim) # A little trick here so I can loop through the data frame generated below for what I want to study
assign("None", NULL)


# Kula vs. None_Kula within Oceania
# Make a dataframe to run the functions

G <- c("Kula","Non_Kula")
E <- c("None","Kula")
n <- length(G)
Hst_d <- tibble(Group=rep(G,7),
                Exclude=rep(E,7),
                Begin=c(rep(1,n),rep(2,n),rep(3,n),rep(4,n),rep(5,n),rep(6,n),rep(10,n)),
                End=c(rep(5,n),rep(6,n),rep(7,n),rep(8,n),rep(9,n),rep(10,n),rep(Inf,n)),
                Hst_L_mean=rep(0,7*n),
                Hst_L_std=rep(0,7*n),
)

# Run functions for each group and length range
for (i in 1:nrow(Hst_d)){
  Begin = Hst_d$Begin[i]
  End = Hst_d$End[i]
  Group = Hst_d$Group[i]
  Exclude = Hst_d$Exclude[i]
  Hst_mean_std_L <- Hst_mean_std(IBD, Info, min=Begin, max=End, all=D, within=eval(as.symbol(Group)), within_exclude=eval(as.symbol(Exclude)))
  Hst_d$Hst_L_mean[i] = Hst_mean_std_L[1]
  Hst_d$Hst_L_std[i] = Hst_mean_std_L[2]
}


# Plot the result
Hst_d$Group <- factor(Hst_d$Group, levels=c("Kula","Non_Kula"), ordered=T)

Hst_d[Hst_d$Begin==10,]$End <- 12

Group_color = c("Kula"="#6A3884", "Non_Kula"="#C4B286")
Hst_d %>% ggplot() +
  geom_point(aes(x = End, y = Hst_L_mean, color=Group), size = 4) +
  geom_errorbar(aes(x = End, y = Hst_L_mean, ymin = Hst_L_mean - 3 * Hst_L_std, ymax = Hst_L_mean + 3 * Hst_L_std, color=Group), width=.5) +
  geom_smooth(aes(x = End, y = Hst_L_mean, color = Group), linetype="dotted",  se =F, size=0.75) +
  scale_x_continuous(breaks = c(5, 6, 7, 8, 9, 10, 12), 
                     labels = c("5","6","7","8","9","10","Over 10"),
                     sec.axis = dup_axis(name = "Lower bound of IBD block length (cM)",
                                         labels = c("1","2","3","4","5","6","10"))) +
  scale_color_manual(values = Group_color) +
  labs(x = "Upper bound of IBD block length (cM)", y = "Relative IBD similarity", color = "Group") +
  #facet_wrap(.~Group, nrow=Nrow) +
  theme_bw() +
  theme(legend.text = element_text(size = 12), legend.title = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14)) +
  #theme(legend.position = c(0.75, 0.1)) +
  theme(legend.background = element_rect(color="black", 
                                         linetype="solid"))




# Do a permutation test here. Generate an expected distribution of Hst from randomly assigned Kula/non-Kula groups.
E_Hst_d <- tibble(Group1=rep(G[1],7),
                  Group2=rep(G[2],7),
                  Begin=c(rep(1,1),rep(2,1),rep(3,1),rep(4,1),rep(5,1),rep(6,1),rep(10,1)),
                  End=c(rep(5,1),rep(6,1),rep(7,1),rep(8,1),rep(9,1),rep(10,1),rep(Inf,1)),
                  Hst_L_diff_mean=rep(0,7),
                  Hst_L_diff_std=rep(0,7),
)


total_l <- list()

# Do it 1000 times
for (i in 1:nrow(E_Hst_d)){
  Begin = E_Hst_d$Begin[i]
  End = E_Hst_d$End[i]
  Group1 = E_Hst_d$Group1[i]
  Group2 = E_Hst_d$Group2[i]
  E_Hst_diff_list <- sapply(1:1000, function(x) Sample_Hst_diff(seed=x, IBD=IBD, Info=Info, Measure="L", 
                                                                min=Begin, max=End, all=D, 
                                                                T_all=eval(as.symbol(Group2)), within=eval(as.symbol(Group1))))
  total_l <- c(total_l, E_Hst_diff_list)
  E_Hst_d$Hst_L_diff_mean[i] = mean(E_Hst_diff_list)
  E_Hst_d$Hst_L_diff_std[i] = std(E_Hst_diff_list)
}

# This is the observed ones from real data that we have done above
O_Hst_d <- Hst_d %>% group_by(Begin) %>% mutate(Hst_L_diff=-(Hst_L_mean-lag(Hst_L_mean))) %>% select(Begin, End, Hst_L_diff) %>% drop_na()

# Save the results
# Modify according to your path
save(Hst_d, E_Hst_d, total_l, O_Hst_d, file="/mnt/scratch/dang/Kula/IBD/conHap/Hst/Hst.All_Kula.vs.None_Kula.Oce.Rdata")
#load("/mnt/scratch/dang/Kula/IBD/conHap/Hst/Hst.All_Kula.vs.None_Kula.Oce.Rdata")


# Plot the permutation results
# Make a tibble of the distribution by length range
d1 <- tibble(
  Hst_L_diff = as.numeric(total_l),
  Begin=as.character(c(rep(1,1000),rep(2,1000),rep(3,1000),rep(4,1000),rep(5,1000),rep(6,1000),rep(10,1000))),
  End=as.character(c(rep(5,1000),rep(6,1000),rep(7,1000),rep(8,1000),rep(9,1000),rep(10,1000),rep(Inf,1000)))
)

# Order the begins
d1$Begin <- factor(d1$Begin, levels=c("1","2","3","4","5","6","10"), ordered=T)

# Plot!
d1 %>% ggplot() +
  stat_boxplot_custom(aes(x=Begin, y=Hst_L_diff), qs = c(.05, .25, 0.5, 0.75, 0.95), fill="grey", width=0.4) +
  geom_point(data=O_Hst_d, aes(x=as.character(Begin), y=Hst_L_diff), fill="red", size=4, pch=23, alpha=0.75) +
  scale_x_discrete(labels=c("1"="2.7", "2"="1.5", "3"="1.1", "4"="0.8", "5"="0.7", "6"="0.6", "10"="0.2")) +
  theme(axis.line.x = element_line(color="black", size = 0.5, linetype = 1),
        axis.line.y = element_line(color="black", size = 0.5, linetype = 1)) +
  labs(x = "kya", y = "Difference of relative IBD similarity") +
  theme_bw() +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 16)) +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14)) +
  theme(axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16))
