
################
### This script compares 3 clades with regards to markov chain model, geographical spread, and zeta diversity.
### Run script1, script3, and script4 before running this script.
################

setwd("C:/Users/garrett.janzen/OneDrive - USDA/Projects/IAV_Env_Eco")
load("IAV_Sources_and_Sinks.RData") # required to run this script

##########################

library("ggplot2")
library("gplots")
library("reshape2")
library("markovchain")
library("diagram")
library("lattice")
library("ade4")
library("vegan")
library("igraph")
library("ecodist")
library("zetadiv")

##########################

data_sub <- data_ss[,c("Date","State","UID_complex"),]
rm(plot1, plot1_poster, plot2, plot3, plot3_late, plot3_nolabels, plot3_poster, plot3_PPPPPP, plot4)
rm(seasonal, sink_table, source_table, list_sinks, list_sources, geodistm)
rm(df, df_mut, df_mut_late, df_mut_PPPPPP, df_mut2)
rm(model, nodes, nodes_m, nodes_t, m)
rm(con, const, constellations, constellations_exclude, cooksD_sink, cooksD_source, data_dmax, data_dmin, dmax, dmin)
rm(i, k, n, r, tempts, persistence_minimum, network_threshold, Maxdist, statestring)

data_H32010 <- data_sub[grep("H3-2010", data_sub$UID_complex), ];dim(data_H32010)
data_H1A333 <- data_sub[grep("H1-1A.3.3.3", data_sub$UID_complex), ];dim(data_H1A333)
data_H1A113 <- data_sub[grep("H1-1A.1.1.3", data_sub$UID_complex), ];dim(data_H1A113)
UIDs_complex_H32010 <- unique(data_H32010$UID_complex)
UIDs_complex_H1A333 <- unique(data_H1A333$UID_complex)
UIDs_complex_H1A113 <- unique(data_H1A113$UID_complex)

mcstatesubset <- c("(IA|MN|NC|IL|IN|NE|MO|OH|OK|SD)")

######################## H32010 ########################

statestringlist <- list()
statestringlist_repeat <- list()
statestring_repeat <- i <- r <- k <- NULL

for(i in 1:length(UIDs_complex_H32010)){
  const <- UIDs_complex_H32010[i]
  dat <- data_H32010[which(data_H32010$UID_complex == const),];dim(dat)
  statestring_repeat <- dat$State
  statestringlist_repeat <- append(statestringlist_repeat, list(statestring_repeat))
  names(statestringlist_repeat)[i] <- const
  statestring_repeat <- NULL
  
  dat <- dat[which(!duplicated(dat$State)),];dim(dat)
  statestring <- dat$State
  statestringlist <- append(statestringlist, list(statestring))
  names(statestringlist)[i] <- const
  statestring <- NULL
}

statestringlist_xx <- statestringlist
statestringlist_repeat_xx <- statestringlist_repeat
data_dmax <- as.Date(max(data$Date))

for(i in 1:length(statestringlist_xx)){
  const <- names(statestringlist_xx)[i]
  
  dat <- data_H32010[which(data_H32010$UID_complex == const),];dim(dat)
  # Below, we append "XX" to the state string if the constellation's final detection is before the cut-off date, with uncertainty window
  xxstring <- ifelse(max(dat$Date) < data_dmax-window, list(c(statestringlist_repeat_xx[[i]], "XX")), list(c(statestringlist_repeat_xx[[i]]))) ### GGG < or > ?
  statestringlist_repeat_xx[[i]] <- unlist(xxstring)
  xxstring <- NULL
  
  dat <- dat[which(!duplicated(dat$State)),];dim(dat)
  xxstring <- ifelse(max(dat$Date) < data_dmax-window, list(c(statestringlist_xx[[i]], "XX")), list(c(statestringlist_xx[[i]]))) ### GGG < or > ?
  statestringlist_xx[[i]] <- unlist(xxstring)
  xxstring <- NULL
}

m <- do.call("rbind", lapply(statestringlist, function(x) cbind(head(x, -1), tail(x, -1))))
mc <- markovchainFit(m)
est <- mc$estimate
tm <- est@transitionMatrix

m_repeat <- do.call("rbind", lapply(statestringlist_repeat, function(x) cbind(head(x, -1), tail(x, -1))))
mc_repeat <- markovchainFit(m_repeat)
est_repeat <- mc_repeat$estimate
tm_repeat_H32010 <- est_repeat@transitionMatrix

m_xx <- do.call("rbind", lapply(statestringlist_xx, function(x) cbind(head(x, -1), tail(x, -1))))
mc_xx <- markovchainFit(m_xx)
est_xx <- mc_xx$estimate
tm_xx <- est_xx@transitionMatrix

m_repeat_xx <- do.call("rbind", lapply(statestringlist_repeat_xx, function(x) cbind(head(x, -1), tail(x, -1))))
mc_repeat_xx <- markovchainFit(m_repeat_xx)
est_repeat_xx <- mc_repeat_xx$estimate
tm_repeat_xx <- est_repeat_xx@transitionMatrix

tm_melt <- melt(tm)
colnames(tm_melt) <- c("Origin", "Destination", "Probability")
tm_melt_repeat_H32010 <- melt(tm_repeat_H32010)
colnames(tm_melt_repeat_H32010) <- c("Origin", "Destination", "Probability")

hmplot_repeat_H32010 <- ggplot(tm_melt_repeat_H32010, aes(Destination, Origin, fill = Probability)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low="#CCFFFF", high="dodgerblue4") +
  geom_tile() +
  labs(title="Markov chain transition matrix");hmplot_repeat_H32010

statestringlist_subset <- list()
statestringlist_repeat_subset <- list()

for(i in 1:length(UIDs_complex_H32010)){
  # i <- 14
  const <- UIDs_complex_H32010[i]
  dat <- data_H32010[which(data_H32010$UID_complex == const),];dim(dat)
  statestring_repeat_subset <- dat$State
  ### To focus on a subset of states, filter here:
  statestring_repeat_subset <- statestring_repeat_subset[grepl(paste0(mcstatesubset), statestring_repeat_subset)]
  statestringlist_repeat_subset <- append(statestringlist_repeat_subset, list(statestring_repeat_subset))
  names(statestringlist_repeat_subset)[i] <- const
  statestring_repeat_subset <- NULL
  
  dat <- dat[which(!duplicated(dat$State)),];dim(dat)
  statestring_subset <- dat$State
  ### To focus on a subset of states, filter here:
  statestring_subset <- statestring_subset[grepl(paste0(mcstatesubset), statestring_subset)]
  statestringlist_subset <- append(statestringlist_subset, list(statestring_subset))
  names(statestringlist_subset)[i] <- const
  statestring_subset <- NULL

}

m_subset <- do.call("rbind", lapply(statestringlist_subset, function(x) cbind(head(x, -1), tail(x, -1))))
mc_subset <- markovchainFit(m_subset)
est_subset <- mc_subset$estimate
tm_subset_H32010 <- est_subset@transitionMatrix
tm_melt_H32010 <- melt(tm_subset_H32010)
colnames(tm_melt_H32010) <- c("Origin", "Destination", "Probability")
hmplot_subset_H32010 <- ggplot(tm_melt_H32010, aes(Destination, Origin, fill = Probability)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low="#CCFFFF", high="dodgerblue4") +
  geom_tile() +
  labs(title="Markov chain transition matrix");hmplot_subset_H32010

m_repeat_subset <- do.call("rbind", lapply(statestringlist_repeat_subset, function(x) cbind(head(x, -1), tail(x, -1))))
mc_repeat_subset <- markovchainFit(m_repeat_subset)
est_repeat_subset <- mc_repeat_subset$estimate
tm_repeat_subset_H32010 <- est_repeat_subset@transitionMatrix
tm_melt_repeat_subset_H32010 <- melt(tm_repeat_subset_H32010)
colnames(tm_melt_repeat_subset_H32010) <- c("Origin", "Destination", "Probability")
hmplot_repeat_subset_H32010 <- ggplot(tm_melt_repeat_subset_H32010, aes(Destination, Origin, fill = Probability)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low="#CCFFFF", high="dodgerblue4") +
  geom_tile() +
  labs(title="Markov chain transition matrix");hmplot_repeat_subset_H32010



######################## H1A333 ########################

statestringlist <- list()
statestringlist_repeat <- list()
statestring_repeat <- i <- r <- k <- NULL

for(i in 1:length(UIDs_complex_H1A333)){
  const <- UIDs_complex_H1A333[i]
  dat <- data_H1A333[which(data_H1A333$UID_complex == const),];dim(dat)
  statestring_repeat <- dat$State
  statestringlist_repeat <- append(statestringlist_repeat, list(statestring_repeat))
  names(statestringlist_repeat)[i] <- const
  statestring_repeat <- NULL
  
  dat <- dat[which(!duplicated(dat$State)),];dim(dat)
  statestring <- dat$State
  statestringlist <- append(statestringlist, list(statestring))
  names(statestringlist)[i] <- const
  statestring <- NULL
}

statestringlist_xx <- statestringlist
statestringlist_repeat_xx <- statestringlist_repeat
data_dmax <- as.Date(max(data$Date))

for(i in 1:length(statestringlist_xx)){
  const <- names(statestringlist_xx)[i]
  
  dat <- data_H1A333[which(data_H1A333$UID_complex == const),];dim(dat)
  # Below, we append "XX" to the state string if the constellation's final detection is before the cut-off date, with uncertainty window
  xxstring <- ifelse(max(dat$Date) < data_dmax-window, list(c(statestringlist_repeat_xx[[i]], "XX")), list(c(statestringlist_repeat_xx[[i]]))) ### GGG < or > ?
  statestringlist_repeat_xx[[i]] <- unlist(xxstring)
  xxstring <- NULL
  
  dat <- dat[which(!duplicated(dat$State)),];dim(dat)
  xxstring <- ifelse(max(dat$Date) < data_dmax-window, list(c(statestringlist_xx[[i]], "XX")), list(c(statestringlist_xx[[i]]))) ### GGG < or > ?
  statestringlist_xx[[i]] <- unlist(xxstring)
  xxstring <- NULL
}

m <- do.call("rbind", lapply(statestringlist, function(x) cbind(head(x, -1), tail(x, -1))))
mc <- markovchainFit(m)
est <- mc$estimate
tm <- est@transitionMatrix

m_repeat <- do.call("rbind", lapply(statestringlist_repeat, function(x) cbind(head(x, -1), tail(x, -1))))
mc_repeat <- markovchainFit(m_repeat)
est_repeat <- mc_repeat$estimate
tm_repeat_H1A333 <- est_repeat@transitionMatrix

m_xx <- do.call("rbind", lapply(statestringlist_xx, function(x) cbind(head(x, -1), tail(x, -1))))
mc_xx <- markovchainFit(m_xx)
est_xx <- mc_xx$estimate
tm_xx <- est_xx@transitionMatrix

m_repeat_xx <- do.call("rbind", lapply(statestringlist_repeat_xx, function(x) cbind(head(x, -1), tail(x, -1))))
mc_repeat_xx <- markovchainFit(m_repeat_xx)
est_repeat_xx <- mc_repeat_xx$estimate
tm_repeat_xx <- est_repeat_xx@transitionMatrix

tm_melt <- melt(tm)
colnames(tm_melt) <- c("Origin", "Destination", "Probability")
tm_melt_repeat_H1A333 <- melt(tm_repeat_H1A333)
colnames(tm_melt_repeat_H1A333) <- c("Origin", "Destination", "Probability")

hmplot_repeat_H1A333 <- ggplot(tm_melt_repeat_H1A333, aes(Destination, Origin, fill = Probability)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low="#CCFFFF", high="dodgerblue4") +
  geom_tile() +
  labs(title="Markov chain transition matrix");hmplot_repeat_H1A333


statestringlist_subset <- list()
statestringlist_repeat_subset <- list()

for(i in 1:length(UIDs_complex_H1A333)){
  # i <- 14
  const <- UIDs_complex_H1A333[i]
  dat <- data_H1A333[which(data_H1A333$UID_complex == const),];dim(dat)
  statestring_repeat_subset <- dat$State
  ### To focus on a subset of states, filter here:
  statestring_repeat_subset <- statestring_repeat_subset[grepl(paste0(mcstatesubset), statestring_repeat_subset)]
  statestringlist_repeat_subset <- append(statestringlist_repeat_subset, list(statestring_repeat_subset))
  names(statestringlist_repeat_subset)[i] <- const
  statestring_repeat_subset <- NULL
  
  dat <- dat[which(!duplicated(dat$State)),];dim(dat)
  statestring_subset <- dat$State
  ### To focus on a subset of states, filter here:
  statestring_subset <- statestring_subset[grepl(paste0(mcstatesubset), statestring_subset)]
  statestringlist_subset <- append(statestringlist_subset, list(statestring_subset))
  names(statestringlist_subset)[i] <- const
  statestring_subset <- NULL
  
}

m_subset <- do.call("rbind", lapply(statestringlist_subset, function(x) cbind(head(x, -1), tail(x, -1))))
mc_subset <- markovchainFit(m_subset)
est_subset <- mc_subset$estimate
tm_subset_H1A333 <- est_subset@transitionMatrix
tm_melt_H1A333 <- melt(tm_subset_H1A333)
colnames(tm_melt_H1A333) <- c("Origin", "Destination", "Probability")
hmplot_subset_H1A333 <- ggplot(tm_melt_H1A333, aes(Destination, Origin, fill = Probability)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low="#CCFFFF", high="dodgerblue4") +
  geom_tile() +
  labs(title="Markov chain transition matrix");hmplot_subset_H1A333

m_repeat_subset <- do.call("rbind", lapply(statestringlist_repeat_subset, function(x) cbind(head(x, -1), tail(x, -1))))
mc_repeat_subset <- markovchainFit(m_repeat_subset)
est_repeat_subset <- mc_repeat_subset$estimate
tm_repeat_subset_H1A333 <- est_repeat_subset@transitionMatrix
tm_melt_repeat_subset_H1A333 <- melt(tm_repeat_subset_H1A333)
colnames(tm_melt_repeat_subset_H1A333) <- c("Origin", "Destination", "Probability")
hmplot_repeat_subset_H1A333 <- ggplot(tm_melt_repeat_subset_H1A333, aes(Destination, Origin, fill = Probability)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low="#CCFFFF", high="dodgerblue4") +
  geom_tile() +
  labs(title="Markov chain transition matrix");hmplot_repeat_subset_H1A333
######################## H1A113 ########################

statestringlist <- list()
statestringlist_repeat <- list()
statestring_repeat <- i <- r <- k <- NULL

for(i in 1:length(UIDs_complex_H1A113)){
  const <- UIDs_complex_H1A113[i]
  dat <- data_H1A113[which(data_H1A113$UID_complex == const),];dim(dat)
  statestring_repeat <- dat$State
  statestringlist_repeat <- append(statestringlist_repeat, list(statestring_repeat))
  names(statestringlist_repeat)[i] <- const
  statestring_repeat <- NULL
  
  dat <- dat[which(!duplicated(dat$State)),];dim(dat)
  statestring <- dat$State
  statestringlist <- append(statestringlist, list(statestring))
  names(statestringlist)[i] <- const
  statestring <- NULL
}

statestringlist_xx <- statestringlist
statestringlist_repeat_xx <- statestringlist_repeat
data_dmax <- as.Date(max(data$Date))

for(i in 1:length(statestringlist_xx)){
  const <- names(statestringlist_xx)[i]
  
  dat <- data_H1A113[which(data_H1A113$UID_complex == const),];dim(dat)
  # Below, we append "XX" to the state string if the constellation's final detection is before the cut-off date, with uncertainty window
  xxstring <- ifelse(max(dat$Date) < data_dmax-window, list(c(statestringlist_repeat_xx[[i]], "XX")), list(c(statestringlist_repeat_xx[[i]]))) ### GGG < or > ?
  statestringlist_repeat_xx[[i]] <- unlist(xxstring)
  xxstring <- NULL
  
  dat <- dat[which(!duplicated(dat$State)),];dim(dat)
  xxstring <- ifelse(max(dat$Date) < data_dmax-window, list(c(statestringlist_xx[[i]], "XX")), list(c(statestringlist_xx[[i]]))) ### GGG < or > ?
  statestringlist_xx[[i]] <- unlist(xxstring)
  xxstring <- NULL
}

m <- do.call("rbind", lapply(statestringlist, function(x) cbind(head(x, -1), tail(x, -1))))
mc <- markovchainFit(m)
est <- mc$estimate
tm <- est@transitionMatrix

m_repeat <- do.call("rbind", lapply(statestringlist_repeat, function(x) cbind(head(x, -1), tail(x, -1))))
mc_repeat <- markovchainFit(m_repeat)
est_repeat <- mc_repeat$estimate
tm_repeat_H1A113 <- est_repeat@transitionMatrix

m_xx <- do.call("rbind", lapply(statestringlist_xx, function(x) cbind(head(x, -1), tail(x, -1))))
mc_xx <- markovchainFit(m_xx)
est_xx <- mc_xx$estimate
tm_xx <- est_xx@transitionMatrix

m_repeat_xx <- do.call("rbind", lapply(statestringlist_repeat_xx, function(x) cbind(head(x, -1), tail(x, -1))))
mc_repeat_xx <- markovchainFit(m_repeat_xx)
est_repeat_xx <- mc_repeat_xx$estimate
tm_repeat_xx <- est_repeat_xx@transitionMatrix

tm_melt <- melt(tm)
colnames(tm_melt) <- c("Origin", "Destination", "Probability")
tm_melt_repeat_H1A113 <- melt(tm_repeat_H1A113)
colnames(tm_melt_repeat_H1A113) <- c("Origin", "Destination", "Probability")

hmplot_repeat_H1A113 <- ggplot(tm_melt_repeat_H1A113, aes(Destination, Origin, fill = Probability)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low="#CCFFFF", high="dodgerblue4") +
  geom_tile() +
  labs(title="Markov chain transition matrix");hmplot_repeat_H1A113

statestringlist_subset <- list()
statestringlist_repeat_subset <- list()

for(i in 1:length(UIDs_complex_H1A113)){
  # i <- 14
  const <- UIDs_complex_H1A113[i]
  dat <- data_H1A113[which(data_H1A113$UID_complex == const),];dim(dat)
  statestring_repeat_subset <- dat$State
  ### To focus on a subset of states, filter here:
  statestring_repeat_subset <- statestring_repeat_subset[grepl(paste0(mcstatesubset), statestring_repeat_subset)]
  statestringlist_repeat_subset <- append(statestringlist_repeat_subset, list(statestring_repeat_subset))
  names(statestringlist_repeat_subset)[i] <- const
  statestring_repeat_subset <- NULL
  
  dat <- dat[which(!duplicated(dat$State)),];dim(dat)
  statestring_subset <- dat$State
  ### To focus on a subset of states, filter here:
  statestring_subset <- statestring_subset[grepl(paste0(mcstatesubset), statestring_subset)]
  statestringlist_subset <- append(statestringlist_subset, list(statestring_subset))
  names(statestringlist_subset)[i] <- const
  statestring_subset <- NULL
  
}

m_subset <- do.call("rbind", lapply(statestringlist_subset, function(x) cbind(head(x, -1), tail(x, -1))))
mc_subset <- markovchainFit(m_subset)
est_subset <- mc_subset$estimate
tm_subset_H1A113 <- est_subset@transitionMatrix
tm_melt_H1A113 <- melt(tm_subset_H1A113)
colnames(tm_melt_H1A113) <- c("Origin", "Destination", "Probability")
hmplot_subset_H1A113 <- ggplot(tm_melt_H1A113, aes(Destination, Origin, fill = Probability)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low="#CCFFFF", high="dodgerblue4") +
  geom_tile() +
  labs(title="Markov chain transition matrix");hmplot_subset_H1A113

m_repeat_subset <- do.call("rbind", lapply(statestringlist_repeat_subset, function(x) cbind(head(x, -1), tail(x, -1))))
mc_repeat_subset <- markovchainFit(m_repeat_subset)
est_repeat_subset <- mc_repeat_subset$estimate
tm_repeat_subset_H1A113 <- est_repeat_subset@transitionMatrix
tm_melt_repeat_subset_H1A113 <- melt(tm_repeat_subset_H1A113)
colnames(tm_melt_repeat_subset_H1A113) <- c("Origin", "Destination", "Probability")
hmplot_repeat_subset_H1A113 <- ggplot(tm_melt_repeat_subset_H1A113, aes(Destination, Origin, fill = Probability)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low="#CCFFFF", high="dodgerblue4") +
  geom_tile() +
  labs(title="Markov chain transition matrix");hmplot_repeat_subset_H1A113

###
hmplot_repeat_H32010
hmplot_repeat_H1A113
hmplot_repeat_H1A333

hmplot_subset_H32010
hmplot_subset_H1A113
hmplot_subset_H1A333

hmplot_repeat_subset_H32010
hmplot_repeat_subset_H1A113
hmplot_repeat_subset_H1A333
###

all_states <- unique(c(rownames(tm_repeat_H32010), rownames(tm_repeat_H1A113), rownames(tm_repeat_H1A333)))
subset_states <- c("IA","MN","NC","IL","IN","NE","MO","OH","OK","SD")

expand_matrix <- function(mat, all_states) {
  mat_expanded <- matrix(0, nrow=length(all_states), ncol=length(all_states), 
                         dimnames = list(all_states, all_states))
  common_states <- intersect(rownames(mat), all_states)
  common_states_col <- intersect(colnames(mat), all_states)
  mat_expanded[common_states, common_states_col] <- mat[common_states, common_states_col]
  return(mat_expanded)
}

tm_repeat_H32010_expanded <- expand_matrix(tm_repeat_H32010, all_states);dim(tm_repeat_H32010_expanded)
tm_repeat_H32010_expanded <- tm_repeat_H32010_expanded[order(rownames(tm_repeat_H32010_expanded)), order(colnames(tm_repeat_H32010_expanded))]
tm_repeat_H32010_expanded_melt <- melt(tm_repeat_H32010_expanded)
colnames(tm_repeat_H32010_expanded_melt) <- c("Origin", "Destination", "Probability")

tm_repeat_H1A113_expanded <- expand_matrix(tm_repeat_H1A113, all_states);dim(tm_repeat_H1A113_expanded)
tm_repeat_H1A113_expanded <- tm_repeat_H1A113_expanded[order(rownames(tm_repeat_H1A113_expanded)), order(colnames(tm_repeat_H1A113_expanded))]
tm_repeat_H1A113_expanded_melt <- melt(tm_repeat_H1A113_expanded)
colnames(tm_repeat_H1A113_expanded_melt) <- c("Origin", "Destination", "Probability")

tm_repeat_H1A333_expanded <- expand_matrix(tm_repeat_H1A333, all_states);dim(tm_repeat_H1A333_expanded)
tm_repeat_H1A333_expanded <- tm_repeat_H1A333_expanded[order(rownames(tm_repeat_H1A333_expanded)), order(colnames(tm_repeat_H1A333_expanded))]
tm_repeat_H1A333_expanded_melt <- melt(tm_repeat_H1A333_expanded)
colnames(tm_repeat_H1A333_expanded_melt) <- c("Origin", "Destination", "Probability")

hmplot_repeat_H32010 <- ggplot(tm_repeat_H32010_expanded_melt, aes(Destination, Origin, fill = Probability)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low="#CCFFFF", high="dodgerblue4") +
  geom_tile() +
  labs(title="Markov chain transition matrix");hmplot_repeat_H32010
hmplot_repeat_H1A113 <- ggplot(tm_repeat_H1A113_expanded_melt, aes(Destination, Origin, fill = Probability)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low="#CCFFFF", high="dodgerblue4") +
  geom_tile() +
  labs(title="Markov chain transition matrix");hmplot_repeat_H1A113
hmplot_repeat_H1A333 <- ggplot(tm_repeat_H1A333_expanded_melt, aes(Destination, Origin, fill = Probability)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low="#CCFFFF", high="dodgerblue4") +
  geom_tile() +
  labs(title="Markov chain transition matrix");hmplot_repeat_H1A333

png("Plots/script5_hmplot_repeat_H32010.png", width=650);hmplot_repeat_H32010;dev.off()
png("Plots/script5_hmplot_repeat_H1A113.png", width=650);hmplot_repeat_H1A113;dev.off()
png("Plots/script5_hmplot_repeat_H1A333.png", width=650);hmplot_repeat_H1A333;dev.off()

tm_subset_H32010_expanded <- expand_matrix(tm_subset_H32010, subset_states);dim(tm_subset_H32010_expanded)
tm_subset_H32010_expanded <- tm_subset_H32010_expanded[order(rownames(tm_subset_H32010_expanded)), order(colnames(tm_subset_H32010_expanded))]
tm_subset_H32010_expanded_melt <- melt(tm_subset_H32010_expanded)
colnames(tm_subset_H32010_expanded_melt) <- c("Origin", "Destination", "Probability")

tm_subset_H1A333_expanded <- expand_matrix(tm_subset_H1A333, subset_states);dim(tm_subset_H1A333_expanded)
tm_subset_H1A333_expanded <- tm_subset_H1A333_expanded[order(rownames(tm_subset_H1A333_expanded)), order(colnames(tm_subset_H1A333_expanded))]
tm_subset_H1A333_expanded_melt <- melt(tm_subset_H1A333_expanded)
colnames(tm_subset_H1A333_expanded_melt) <- c("Origin", "Destination", "Probability")

tm_subset_H1A113_expanded <- expand_matrix(tm_subset_H1A113, subset_states);dim(tm_subset_H1A113_expanded)
tm_subset_H1A113_expanded <- tm_subset_H1A113_expanded[order(rownames(tm_subset_H1A113_expanded)), order(colnames(tm_subset_H1A113_expanded))]
tm_subset_H1A113_expanded_melt <- melt(tm_subset_H1A113_expanded)
colnames(tm_subset_H1A113_expanded_melt) <- c("Origin", "Destination", "Probability")

hmplot_subset_H32010 <- ggplot(tm_subset_H32010_expanded_melt, aes(Destination, Origin, fill = Probability)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low="#CCFFFF", high="dodgerblue4") +
  geom_tile() +
  labs(title="Markov chain transition matrix");hmplot_subset_H32010
hmplot_subset_H1A113 <- ggplot(tm_subset_H1A113_expanded_melt, aes(Destination, Origin, fill = Probability)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low="#CCFFFF", high="dodgerblue4") +
  geom_tile() +
  labs(title="Markov chain transition matrix");hmplot_subset_H1A113
hmplot_subset_H1A333 <- ggplot(tm_subset_H1A333_expanded_melt, aes(Destination, Origin, fill = Probability)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low="#CCFFFF", high="dodgerblue4") +
  geom_tile() +
  labs(title="Markov chain transition matrix");hmplot_subset_H1A333

tm_repeat_subset_H32010_expanded <- expand_matrix(tm_repeat_subset_H32010, subset_states);dim(tm_repeat_subset_H32010_expanded)
tm_repeat_subset_H32010_expanded <- tm_repeat_subset_H32010_expanded[order(rownames(tm_repeat_subset_H32010_expanded)), order(colnames(tm_repeat_subset_H32010_expanded))]
tm_repeat_subset_H32010_expanded_melt <- melt(tm_repeat_subset_H32010_expanded)
colnames(tm_repeat_subset_H32010_expanded_melt) <- c("Origin", "Destination", "Probability")

tm_repeat_subset_H1A333_expanded <- expand_matrix(tm_repeat_subset_H1A333, subset_states);dim(tm_repeat_subset_H1A333_expanded)
tm_repeat_subset_H1A333_expanded <- tm_repeat_subset_H1A333_expanded[order(rownames(tm_repeat_subset_H1A333_expanded)), order(colnames(tm_repeat_subset_H1A333_expanded))]
tm_repeat_subset_H1A333_expanded_melt <- melt(tm_repeat_subset_H1A333_expanded)
colnames(tm_repeat_subset_H1A333_expanded_melt) <- c("Origin", "Destination", "Probability")

tm_repeat_subset_H1A113_expanded <- expand_matrix(tm_repeat_subset_H1A113, subset_states);dim(tm_repeat_subset_H1A113_expanded)
tm_repeat_subset_H1A113_expanded <- tm_repeat_subset_H1A113_expanded[order(rownames(tm_repeat_subset_H1A113_expanded)), order(colnames(tm_repeat_subset_H1A113_expanded))]
tm_repeat_subset_H1A113_expanded_melt <- melt(tm_repeat_subset_H1A113_expanded)
colnames(tm_repeat_subset_H1A113_expanded_melt) <- c("Origin", "Destination", "Probability")

hmplot_repeat_subset_H32010 <- ggplot(tm_repeat_subset_H32010_expanded_melt, aes(Destination, Origin, fill = Probability)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low="#CCFFFF", high="dodgerblue4") +
  geom_tile() +
  labs(title="Markov chain transition matrix");hmplot_repeat_subset_H32010
hmplot_repeat_subset_H1A113 <- ggplot(tm_repeat_subset_H1A113_expanded_melt, aes(Destination, Origin, fill = Probability)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low="#CCFFFF", high="dodgerblue4") +
  geom_tile() +
  labs(title="Markov chain transition matrix");hmplot_repeat_subset_H1A113
hmplot_repeat_subset_H1A333 <- ggplot(tm_repeat_subset_H1A333_expanded_melt, aes(Destination, Origin, fill = Probability)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low="#CCFFFF", high="dodgerblue4") +
  geom_tile() +
  labs(title="Markov chain transition matrix");hmplot_repeat_subset_H1A333

##########################

m_control <- vegan::mantel(as.dist(tm_repeat_subset_H32010_expanded),as.dist(tm_repeat_subset_H32010_expanded));m_control
set.seed(0)
m_H32010_H1A333 <- vegan::mantel(as.dist(tm_repeat_subset_H32010_expanded),
                                as.dist(tm_repeat_subset_H1A333_expanded));m_H32010_H1A333
set.seed(0)
m_H32010_H1A113 <- vegan::mantel(as.dist(tm_repeat_subset_H32010_expanded),
                                as.dist(tm_repeat_subset_H1A113_expanded));m_H32010_H1A113
set.seed(0)
m_H1A113_H1A333 <- vegan::mantel(as.dist(tm_repeat_subset_H1A113_expanded),
                                as.dist(tm_repeat_subset_H1A333_expanded));m_H1A113_H1A333

###########################
###########################
###########################

# https://datascience.blog.wzb.eu/2018/05/31/three-ways-of-visualizing-a-graph-on-a-map/

colnames(centroids)[2] <- "lon"
centroids <- merge(centroids, data_usda_agg, all.x=TRUE, all.y=FALSE, by.x="postal_code", by.y="State")
centroids$state <- tolower(centroids$state)
centroids_formerge <- centroids[,4:5]
colnames(centroids_formerge) <- c("state","Swine Inventory (million)")
centroids_formerge$`Swine Inventory (million)` <- centroids_formerge$`Swine Inventory (million)`/1000000

tm_list <- list(tm_subset_H32010_expanded,tm_subset_H1A333_expanded,tm_subset_H1A113_expanded)
name_list <- list("H32010","H1A333","H1A113")
hist(tm_subset_H32010_expanded)
hist(tm_subset_H1A333_expanded)
hist(tm_subset_H1A113_expanded)
crit_thresh <- 0.30

for(i in 1:length(tm_list)){
  tm_map <- tm_list[[i]]
  name <- name_list[[i]]
  tm_try <- as.data.frame(c(tm_map))
  colnames(tm_try) <- "probability"
  tm_try$to <- rep(rownames(tm_map),each=ncol(tm_map))
  tm_try$from <- rep(colnames(tm_map),nrow(tm_map))
  tm_try2 <- merge(tm_try, centroids, all.x=TRUE, all.y=FALSE, by.x="from", by.y="postal_code");head(tm_try2)
  tm_try2$Production <- NULL
  colnames(tm_try2)[4:5] <- c("y","x")
  tm_try3 <- merge(tm_try2, centroids, all.x=TRUE, all.y=FALSE, by.x="to", by.y="postal_code");head(tm_try3)
  tm_try3$Production <- NULL
  colnames(tm_try3)[7:8] <- c("yend","xend")
  tm_try4 <- tm_try3[tm_try3$y!=tm_try3$yend,]
  tm_try4 <- tm_try4[tm_try4$x!=tm_try4$xend,]
  # tm_plot <- plot(tm_try4$probability, size=2);tm_plot
  tm_try5 <- tm_try4[tm_try4$probability >= crit_thresh,] #if the probability is below the critical threshold defined above, drop it.
  tm_try5$Production.y <- NULL
  rm(tm_try, tm_try2, tm_try3, tm_try4)
  
  tm_try5$Category <- "other"
  tm_try5$Category <- ifelse(tm_try5$to == "IA", "to Iowa", tm_try5$Category)
  tm_try5$Category <- ifelse(tm_try5$from == "IA", "from Iowa", tm_try5$Category)
  tm_try5$Category <- factor(tm_try5$Category, levels=c("to Iowa","from Iowa","other"))
  colnames(tm_try5)[3] <- "Probability"
  
  centroids_try <- centroids[which(centroids$postal_code %in% tm_try5$to | centroids$postal_code %in% tm_try5$from),]
  g <- graph_from_data_frame(tm_try5, directed = FALSE, vertices = centroids_try)
  
  all_states  <- map_data('state')
  all_states$id  <- 1:nrow(all_states)
  stateData <- merge(all_states,centroids_formerge,by.x="region",by.y="state",all.x=TRUE, all.y=FALSE)
  stateData <- stateData[order(stateData$id), ]
  table(stateData$`Swine Inventory (million)`)#;plot(stateData$`Swine Inventory (million)`)
  
  usa_shapes <- geom_polygon(data = stateData, aes(x = long, y = lat, group = group, fill = `Swine Inventory (million)`), size = 0.15)
  centroids_try$Production = degree(g)
  
  mapcoords <- coord_fixed(xlim = c(-103, -76.5), ylim = c(33, 48))
  
  maptheme <- theme(panel.grid = element_blank()) + 
    theme(axis.text = element_blank()) +
    theme(axis.ticks = element_blank()) +
    theme(axis.title = element_blank()) +
    theme(legend.position = "bottom") +
    theme(panel.grid = element_blank()) +
    theme(panel.background = element_rect(fill = "#596673")) +
    theme(plot.margin = unit(c(0, 0, 0.5, 0), 'cm')) +
    theme(legend.title = element_text(size = 17), 
          legend.text = element_text(size = 10))
  
  pdf(paste0("Plots/script5_map1_probability_", name ,".pdf"), height=7.2, width=12)
  print(ggplot(centroids_try) + usa_shapes + maptheme + mapcoords + 
    scale_fill_gradient(name="Swine inventory (million)",low="grey30",high="pink3", breaks=c(0,3,6,9,12)) +
    geom_curve(data=tm_try5, arrow = arrow(length = unit(0.03, "npc"), type="closed"),
               aes(x = x, y = y, xend = xend, yend = yend, size=Probability, color=Probability),
               curvature = 0.33,
               alpha = 0.85, show.legend=TRUE) +
    scale_color_gradient(name="Transition probability",low="#CCFFFF", high="dodgerblue4") +
    scale_size_continuous(guide = FALSE, range = c(1, 4)) + 
    geom_text(data = centroids_try, aes(x = lon, y = lat, label = postal_code), 
              hjust = 0, nudge_x = 0.2, nudge_y = -0.2,
              size = 7, color = "black", fontface = "bold"))
  dev.off()
  
  png(paste0("Plots/script5_map1_probability_", name ,".png"), height=720, width=1200)
  print(ggplot(centroids_try) + usa_shapes + maptheme + mapcoords + 
          scale_fill_gradient(name="Swine inventory (million)",low="grey30",high="pink3", breaks=c(0,3,6,9,12)) +
          geom_curve(data=tm_try5, arrow = arrow(length = unit(0.03, "npc"), type="closed"),
                     aes(x = x, y = y, xend = xend, yend = yend, size=Probability, color=Probability),
                     curvature = 0.33,
                     alpha = 0.85, show.legend=TRUE) +
          scale_color_gradient(name="Transition probability",low="#CCFFFF", high="dodgerblue4") +
          scale_size_continuous(guide = FALSE, range = c(1, 4)) + 
          geom_text(data = centroids_try, aes(x = lon, y = lat, label = postal_code), 
                    hjust = 0, nudge_x = 0.2, nudge_y = -0.2,
                    size = 7, color = "black", fontface = "bold"))
  dev.off()

  scm_direction <- c("#C8102E", "gold", "gray75")
  scm_direction <- ifelse(length(unique(tm_try5$Category)) == 2 & unique(tm_try5$Category) == c("to Iowa","other"),
                          c("#C8102E","gray75"), scm_direction)

  pdf(paste0("Plots/script5_map2_direction_", name ,".pdf"), height=7.2, width=12)
  print(ggplot(centroids_try) + usa_shapes + maptheme + mapcoords +
    scale_fill_gradient(name="Swine inventory (million)",low="grey30",high="pink3") +
    geom_curve(data=tm_try5, arrow = arrow(length = unit(0.03, "npc"), type="closed"),
               aes(x = x, y = y, xend = xend, yend = yend, size=Probability, color=Category),
               curvature = 0.33,
               alpha = 0.85, show.legend=TRUE) +
    scale_color_manual(values=scm_direction) +
    guides(color = guide_legend(override.aes = list(size = 20))) +
    scale_size_continuous(guide = FALSE, range = c(1, 4)) + 
    geom_text(data = centroids_try, aes(x = lon, y = lat, label = postal_code), 
              hjust = 0, nudge_x = 0.2, nudge_y = -0.2,
              size = 7, color = "black", fontface = "bold"))
  dev.off()
  
  png(paste0("Plots/script5_map2_direction_", name ,".png"), height=720, width=1000)
  print(ggplot(centroids_try) + usa_shapes + maptheme + mapcoords +
          scale_fill_gradient(name="Swine inventory (million)",low="grey30",high="pink3") +
          geom_curve(data=tm_try5, arrow = arrow(length = unit(0.03, "npc"), type="closed"),
                     aes(x = x, y = y, xend = xend, yend = yend, size=Probability, color=Category),
                     curvature = 0.33,
                     alpha = 0.85, show.legend=TRUE) +
          scale_color_manual(values=scm_direction) +
          guides(color = guide_legend(override.aes = list(size = 20))) +
          scale_size_continuous(guide = FALSE, range = c(1, 4)) + 
          geom_text(data = centroids_try, aes(x = lon, y = lat, label = postal_code), 
                    hjust = 0, nudge_x = 0.2, nudge_y = -0.2,
                    size = 7, color = "black", fontface = "bold"))
  dev.off()
}

##########################
# Read in the data, make appropriate data objects
data_H32010 <- data_sub[grep("H3-2010", data_sub$UID_complex), ];dim(data_H32010)
data_H1A333 <- data_sub[grep("H1-1A.3.3.3", data_sub$UID_complex), ];dim(data_H1A333)
data_H1A113 <- data_sub[grep("H1-1A.1.1.3", data_sub$UID_complex), ];dim(data_H1A113)

#########################

onehot_fun <- function(x){
  ifelse(x != 0, 1, 0)}

table_H32010 <- as.data.frame.matrix(table(data_H32010[,c("State","UID_complex")]))
onehot_H32010 <- t(data.frame(lapply(table_H32010,onehot_fun)));colnames(onehot_H32010) <- row.names(table_H32010)
table_H1A333 <- as.data.frame.matrix(table(data_H1A333[,c("State","UID_complex")]))
onehot_H1A333 <- t(data.frame(lapply(table_H1A333,onehot_fun)));colnames(onehot_H1A333) <- row.names(table_H1A333)
table_H1A113 <- as.data.frame.matrix(table(data_H1A113[,c("State","UID_complex")]))
onehot_H1A113 <- t(data.frame(lapply(table_H1A113,onehot_fun)));colnames(onehot_H1A113) <- row.names(table_H1A113)

onehot_H32010_t <- as.data.frame(t(onehot_H32010));onehot_H32010_t$postal_code <- row.names(onehot_H32010_t)
onehot_H32010_t <- merge(centroids[,1:3],onehot_H32010_t, by="postal_code", all.y=TRUE, all.x=FALSE)
onehot_H1A333_t <- as.data.frame(t(onehot_H1A333));onehot_H1A333_t$postal_code <- row.names(onehot_H1A333_t)
onehot_H1A333_t <- merge(centroids[,1:3],onehot_H1A333_t, by="postal_code", all.y=TRUE, all.x=FALSE)
onehot_H1A113_t <- as.data.frame(t(onehot_H1A113));onehot_H1A113_t$postal_code <- row.names(onehot_H1A113_t)
onehot_H1A113_t <- merge(centroids[,1:3],onehot_H1A113_t, by="postal_code", all.y=TRUE, all.x=FALSE)

zeta.decline_H32010 <- Zeta.decline.mc(onehot_H32010_t[,-c(1:3)], onehot_H32010_t[,2:3], orders = 1:nrow(onehot_H32010_t), sam = 100, NON = TRUE)
zeta.decline_H1A333 <- Zeta.decline.mc(onehot_H1A333_t[,-c(1:3)], onehot_H1A333_t[,2:3], orders = 1:nrow(onehot_H1A333_t), sam = 100, NON = TRUE)
zeta.decline_H1A113 <- Zeta.decline.mc(onehot_H1A113_t[,-c(1:3)], onehot_H1A113_t[,2:3], orders = 1:nrow(onehot_H1A113_t), sam = 100, NON = TRUE)

pdf("Plots/script5_zeta_decline_H32010.pdf", height=3, width=9);Plot.zeta.decline(zeta.decline_H32010);dev.off()
pdf("Plots/script5_zeta_decline_H1A333.pdf", height=3, width=9);Plot.zeta.decline(zeta.decline_H1A333);dev.off()
pdf("Plots/script5_zeta_decline_H1A113.pdf", height=3, width=9);Plot.zeta.decline(zeta.decline_H1A113);dev.off()

####
val <- max(c(nrow(table_H32010),nrow(table_H1A333),nrow(table_H1A113)));val
length(zeta.decline_H32010$zeta.val) <- length(zeta.decline_H32010$zeta.val.sd) <- val
length(zeta.decline_H1A333$zeta.val) <- length(zeta.decline_H1A333$zeta.val.sd) <- val
length(zeta.decline_H1A113$zeta.val) <- length(zeta.decline_H1A113$zeta.val.sd) <- val

test <- cbind(
  c(1:val), 
  zeta.decline_H32010$zeta.val,
  zeta.decline_H1A333$zeta.val,
  zeta.decline_H1A113$zeta.val)
colnames(test) <- c("order","H32010","H1A333","H1A113")
test_melt <- melt(test)
colnames(test_melt) <- c("Zeta_Order","Clade","Zeta")

zeta.decline_H32010$aic #exp
zeta.decline_H1A333$aic #exp
zeta.decline_H1A113$aic #exp

test_sd <- cbind(
  zeta.decline_H32010$zeta.order, 
  zeta.decline_H32010$zeta.val.sd,
  zeta.decline_H1A333$zeta.val.sd,
  zeta.decline_H1A113$zeta.val.sd)
colnames(test_sd) <- c("order","H32010","H1A333","H1A113")
test_melt_sd <- melt(test_sd)
colnames(test_melt_sd) <- c("Zeta_Order","Clade","SD")

test_merge_full <- merge(test_melt, test_melt_sd, by=c("Zeta_Order","Clade"))
test_merge_full <- test_merge_full[which(test_merge_full$Clade != "order"),]
test_merge <- test_merge_full[which(test_merge_full$Clade != "full"),]

plot2_full <- ggplot(test_merge_full, aes(y=Zeta, x=Zeta_Order, group=Clade, color=Clade)) +
  geom_line(aes(color=Clade), size=1.5) + 
  geom_point(aes(color=Clade), size=2) +
  xlim(0.95,16) +
  labs(title="",
       x ="Zeta Order",
       y = "Zeta Diversity") +
  # geom_errorbar(aes(ymin=Zeta-SD, ymax=Zeta+SD), size=0.5, width=.3, position=position_dodge(0.1)) + #un-comment for SD bars
  theme(legend.title=element_blank(),
        axis.text.y = element_text(size = 22),
        axis.text.x = element_text(size = 22),
        text = element_text(size = 26));plot2_full
