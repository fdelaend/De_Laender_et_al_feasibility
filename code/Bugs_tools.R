
# SET PARAMETERS OF ANALYSIS -------------------------------------------------
cbPalette         <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                   "#0072B2", "#D55E00", "#CC79A7") #colourblind friendly palette
precision         <- 4 #spatial precision with which to match env. and bio data (decimals of lat/long)
period_lenght     <- 5 #length of a period in years
fraction_of_sites <- 0.10 #only genera that are found back in at least this fraction of sites are considered
cutoff_frac       <- 0.05 #only genera that have a local rel. frequency > cutoff_frac*100 %

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

# PREPPING THE CHEM DATA ---------------------------------------------------------------

## 2004: 2000, 2001, 2002, 2003, 2004 ----

sites_2004 <- read.csv("../data/sites2004.csv") %>% # to get locations of the data
  mutate(lat=round(LAT_DD,precision)) %>% #do the rounding
  mutate(long=round(LON_DD,precision)) %>% 
  group_by(SITE_ID) %>% #site's average location
  summarize(mean(lat), mean(long)) %>%
  rename(lat='mean(lat)') %>% 
  rename(long='mean(long)') 

chem_2004  <- read.csv("../data/chem2004.csv") %>%# to get chemistry data
  left_join(sites_2004, by="SITE_ID") %>%
  select(-ends_with("F")) %>%#kill all cols with F at the end
  select(-c("SITE_ID", "DATE_COL", "VISIT_NO", "TEAM_ID", "DATECHEM", "SAMPLED", #kick out nonrelevant info
            "SAMP_LOC", "H", "COLOR", "DAY_SHIP", "COM_LAB", 
            "COM_FLD", "COM_IM", "DIC", "SE", "NA.", "ZN",
            "OH", "CO3", "ALKCALC", "CATSUM", "ANSUM", "SOBC", 
            "IONSTR", "BALANCE", "ORGION", "CONCAL", "CONDHO")) %>%
  rename(calcium=CA) %>%
  rename(magnesium=MG) %>%
  rename(potassium=K) %>%
  rename(ph=PHSTVL) %>%
  rename(ammonia_n=NH4) %>%
  rename(sulfate=SO4) %>%
  rename(nitrate_n=NO3) %>%
  rename(chloride=CL) %>%
  rename(silica=SIO2)%>%
  rename_with(tolower) 

## 2008-2009----

sites_2008_2009 <- read.csv("../data/sites2008-2009.csv") %>% # to get locations of the data
  mutate(lat=round(LAT_DD83,precision)) %>% #do the rounding
  mutate(long=round(LON_DD83,precision)) %>% 
  group_by(SITE_ID) %>% #site's average location
  summarize(mean(lat), mean(long)) %>%
  rename(lat='mean(lat)') %>% 
  rename(long='mean(long)') 

chem_2008_2009  <- read.csv("../data/chem2008-2009.csv") %>%# to get chemistry data
  left_join(sites_2008_2009, by="SITE_ID") %>%
  select(-ends_with("_ALERT")) %>%#kill all cols with F at the end
  select(-c("PUBLICATION_DATE", "UID", "SITE_ID", "DATE_COL", "VISIT_NO", #kick out nonrelevant info
            "SAMPLE_TYPE", "SAM_CODE", "SAMPLE_CAT", "CHEM_COM", "COLOR"))%>%
  rename(calcium=CA) %>%
  rename(magnesium=MG) %>%
  rename(potassium=K) %>%
  rename(ph=PHLAB) %>%
  rename(ammonia_n=NH4) %>%
  rename(sulfate=SO4) %>%
  rename(nitrate_n=NO3) %>%
  rename(chloride=CL) %>%
  rename(silica=SIO2)%>%
  rename_with(tolower) 

## 2013-2014----

sites_2013_2014 <- read.csv("../data/sites2013-2014.csv") %>% # to get locations of the data
  mutate(lat=round(LAT_DD83,precision)) %>% #do the rounding
  mutate(long=round(LON_DD83,precision)) %>% 
  group_by(SITE_ID) %>% #site's average location
  summarize(mean(lat), mean(long)) %>%
  rename(lat='mean(lat)') %>% 
  rename(long='mean(long)') 

chem_2013_2014  <- read.csv("../data/chem2013-2014.csv") %>%# to get chemistry data
  left_join(sites_2013_2014, by="SITE_ID") %>%
  separate(DATE_COL, sep="/", into = c("month", "day", "year")) %>%
  mutate(year=as.numeric(year)) %>%
  select(-contains(c("_BATCH", "_CODE", "_FLAG", "_UNITS", "month", "day",
                     "_ID", "_DATE", "_DD83", "INDEX_", "_TYPE", 
                     "LAB", "UID", "SITETYPE", "DATE_COL", 
                     "VISIT_NO", "AG_ECO9_NM", "MATRIX", "_TIME", 
                     "_MDL", "_RL", "_LRL", "COLOR", "SAMPLE_CAT"))) %>% #kill several cols 
  rename_with(~ gsub("_RESULT", "", .x, fixed = TRUE))%>% #rip off 'results' part
  mutate(across(everything(), ~ as.numeric(.x))) %>%
  rename_with(tolower) 

## 2018-2019 ----

sites_2018_2019 <- read.csv("../data/sites2018-2019.csv") %>% # to get locations of the data
  mutate(lat=round(LAT_DD83,precision)) %>% #do the rounding
  mutate(long=round(LON_DD83,precision)) %>% 
  group_by(SITE_ID) %>% #site's average location
  summarize(mean(lat), mean(long)) %>%
  rename(lat='mean(lat)') %>% 
  rename(long='mean(long)') 

chem_2018_2019  <- read.csv("../data/chem2018-2019.csv") %>%# to get chemistry data
  left_join(sites_2018_2019, by="SITE_ID") %>%
  separate(DATE_COL, sep="/", into = c("month", "day", "year")) %>%
  mutate(year=as.numeric(year)) %>%
  select(-contains(c("_BATCH", "_CODE", "_FLAG", "_UNITS", "month", "day",
                     "_ID", "_DATE", "_DD83", "INDEX_", "_TYPE", 
                     "LAB", "UID", "SITETYPE", "DATE_COL", 
                     "VISIT_NO", "AG_ECO9_NM", "MATRIX", "_TIME", 
                     "_MDL", "_RL", "_LRL", "_COLLECTED", "_FACTOR",
                     "EPA_REG", "AG_ECO9", "STATE", "_VOLUME",
                     "_DISS_", "COLOR", "TKN", "CHLA", "PHEO"))) %>% #kill several cols 
  rename_with(~ gsub("_RESULT", "", .x, fixed = TRUE)) %>% #rip off 'results' part
  mutate(across(everything(), ~ as.numeric(.x))) %>%
  rename_with(tolower)

## COLLAPSE all chemistry data ----------------------------------
chem <- bind_rows(chem_2004, chem_2008_2009, 
                  chem_2013_2014, chem_2018_2019) %>% 
  mutate(period=paste((ceiling(year/period_lenght)-1) * period_lenght,"-",
                      ceiling(year/period_lenght) * period_lenght)) %>% #convert to period
  rename_with(~ paste("chem_", .x, sep = ""), .cols = !any_of(c("year", "lat", "long", "period"))) #give prefix to all chemistru variables for easy manpulation later on

# GET EMPIRICAL INTERACTIONS FROM RGLOBI -------------------------------------------------
# Only those that have observed interactions in rGlobi
invert_key <- read.csv("../data/GeneraTaxonomicInformation.csv") #key to transform  genus  into coarser level

# Only run if trophic_positions.RData not yet available
#library(rglobi)
#trophic_positions   <- NULL #get trophic positions for all genera
#for (i in 1:nrow(invert_key)){
#  trophic_positions <- c(trophic_positions, get_trophic(name=invert_key$Genus[i]))
#  if (i %in% c(1:5)*100) {Sys.sleep(5)}
#}
#save(trophic_positions, file="trophic_positions.RData")
load(file="../data/trophic_positions.RData")
names(trophic_positions) <- invert_key$Genus
all_genera               <- names(trophic_positions)

# LOAD THE INVERT DATA and select the appropriate columns -------------------
invert <- read.csv("../data/invert_dat_density.csv") %>%
  select(all_of(c("Agency", "CollectionYear", "Latitude_dd", "Longitude_dd", "NLCDGroup", all_genera))) %>%
  mutate(stress = case_when(NLCDGroup %in% c("ForestWet", "Other") ~ "low", #, Dichotomy
                            NLCDGroup %in% c("Ag", "Urban") ~ "high")) %>% #, 
  filter(Agency == "EPA") %>% # Only epa data (no env. data at usgs sites)
  select(-Agency) %>% # and remove altogether
  mutate(period=paste((ceiling(CollectionYear/period_lenght)-1) * period_lenght,"-",
                      ceiling(CollectionYear/period_lenght) * period_lenght)) %>% #convert to period
  mutate(lat = round(Latitude_dd, precision)) %>% #round coordinates for matching w env.data
  mutate(long = round(Longitude_dd, precision))
