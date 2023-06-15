
cbPalette         <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                   "#0072B2", "#D55E00", "#CC79A7") #colourblind friendly palette
# FUNCTIONS -------------------------------------------------
# get a span of the data in x (max - min)
span <- function(x) {
  max_min <- as.vector(max(x, na.rm=T) - min(x, na.rm=T))
  ifelse(is.infinite(max_min), NA, max_min)
}

# check if the genera with which a focal genus is interacting are inverts
# data = the output of get_prey_of or get_predators_of (rglobi package); 
# data can contain information on more than one focal genus. In that case, 
# we do the check per focal genus, and only if 
# for at least one of the focal genera the interacting species are invertebrates 
# do we conclude this combined set of genera interacts with invertebrates
check_invert<- function(data) {
  data <- data %>%
    #get out genus per row
    mutate(genus=pmap(., function(source_taxon_name,...) str_split(source_taxon_name,
                                                                   pattern=" ")[[1]][1])) %>%
    unnest(genus) %>%
    #check if interaction is with invert (1 if yes, 0 if not), for every row
    mutate(target_taxon_path_invert = pmap(., function(target_taxon_path, ...) 1*str_detect(target_taxon_path, regex("Arthropoda")))) %>%
    unnest(target_taxon_path_invert) %>%
    #for every genus, compute the average for this score of "invertebrateness"
    group_by(genus) %>%
    summarise(mean_inv = mean(target_taxon_path_invert, na.rm=T))
    return(mean(data$mean_inv, na.rm=T)>0) #only if for at least some of the genera the interacting species are invertebrates do we conclude this combined set of genera interacts with invertebrates 
}

# get the TL of a genus with Globi; 
# name = one or multiple genera, separated by sep. If multiple genera given, it computes the 
# sep = the separator (where the \\ are escapes; otherwise the dot won't work)
get_trophic <- function(name, sep = "\\.", ...){
  genera                 <- str_split(name, sep)[[1]] #if multiple focal genera, split them
  genera_eat_inverts     <- check_invert(get_prey_of(genera)) #does it eat other inverts?
  genera_eaten_by_inverts<- check_invert(get_predators_of(genera)) #does it get eaten by other inverts?
  role                   <- NA #default result, no info on the focal genera
  if ((!is.na(genera_eat_inverts))*(!is.na(genera_eaten_by_inverts))>0){
    if (genera_eat_inverts*(1-genera_eaten_by_inverts)) {role <- "predator"}
    if (genera_eat_inverts*genera_eaten_by_inverts) {role <- "consumer"} 
    if ((1-genera_eat_inverts)*(genera_eaten_by_inverts)) {role <- "resource"}
  }
  return(role)
}

#transform a vector n with abundances of genera 
#to a vector of zeros and ones, where a '1' means the genus is present
#and a '0' if it's absent (aka rel. abundance < cutoff). 
#Represents a unique code for the given community
#If there are NAs in n, these are passed on.
transform_to_taxonomic <- function(n, cutoff=0.05,...){
  return(list(community=paste((n/sum(n, na.rm=T)>cutoff)*1, collapse="_"),
         community_size=sum((n/sum(n)>cutoff_frac))))
}

compute_range <- function(data, ...){
  range <- as.numeric(apply(data, 2, span))
  names(range) <- colnames(data)
  return(list(range=range, 
              n=sum(!is.na(rowSums(data))))) 
}

# Computes the range for each of a series of variables (range_vars),
# observed across various combinations of other variables (combo_vars). 
# The result is pivoted to get the range for every level of 
# the variable pivot_var. These columns containing the computed ranges get the 
# name range_name
compute_range_per_combo <- function(data, combo_vars=c("period", "stress"), 
                                    range_vars, pivot_var="stress", range_name) {
  data %>%
    ungroup() %>%
    select(all_of(c(combo_vars, range_vars))) %>% # Only vars needed 
    nest_by(across({{combo_vars}})) %>%
    ungroup() %>%
    mutate(range = pmap(.,compute_range)) %>%
    unnest_wider(range) %>%
    unnest_wider(range) %>%
    filter(n>2) %>% #kick out combos with less than two observations (none, but just to be complete)
    #select(-all_of(c("data","n"))) %>% #ditch raw data, not needed anymore
    select(-all_of(c("data"))) %>% #ditch raw data, not needed anymore
    pivot_longer(all_of(c(range_vars)), names_to="parameter", values_to = range_name)
}
