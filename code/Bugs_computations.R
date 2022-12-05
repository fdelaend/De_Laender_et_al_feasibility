
# JOIN CHEMISTRY AND INVERTEBRATE DATA -----------------------------
joined_data <- invert %>% 
  left_join(chem, by = c("lat", "long", "period")) %>% #join with chemistry data, based on location and period
  group_by(lat, long, period, stress) %>% #for every combo...
  summarise(across(starts_with("chem_"), ~ mean(.x, na.rm=T)), #... compute the mean of the chemistry parameter
            across(all_of(all_genera), ~ sum(.x, na.rm = TRUE))) %>% #... and the sum of the biology
  ungroup() %>%
  select(-where(~sum(is.na(.x))>3000)) %>% #drop variables with too many NAs or NaNs
  filter(!is.na(stress)) #ditch NA stress levels
  
# SELECT CHEMISTRY VARIABLES THAT DON'T CORRELATE AND REGIONALLY NONRARE TAXA ------------
## identify chemistry variables: run a correlation analysis of the chemistry data to identify correlated variables (more than 0.5) -------
cor_table  <- joined_data %>%
  select(starts_with("chem_")) %>%
  correlate()
cor_table[upper.tri(cor_table)] <- 0
cor_table <- cor_table %>%
  mutate(across(starts_with("chem_"), ~ (abs(.x)<0.5)*1)) %>% 
  rowwise() %>%
  mutate(none_too_high = prod(across(starts_with("chem_")), na.rm = T)) %>% #flag rows with at least one too high correlation
  filter(none_too_high==1)
chemistry_to_keep <- cor_table$term #those that we will use

## identify genera -----------------------
# that have been shown to interact and are not too rare
across_data_set_per_genus    <- colSums(joined_data[,all_genera]>0)
not_too_rare                 <- (across_data_set_per_genus / nrow(joined_data)) > fraction_of_sites
interacting_genera           <- !is.na(trophic_positions)
genera_to_keep               <- all_genera[which(interacting_genera*not_too_rare>0)]

## perform the selection of chemistry and genera to keep,  ----------------
# ...and construct a code for every local community, 
# ...plus ditch the 95-00 data because communities only observed at one site (max two) 
# ...across high stress sites
joined_data_selection <- joined_data %>%
  select(all_of(c("lat", "long", "period", "stress", chemistry_to_keep, genera_to_keep))) %>%
  nest(n=all_of(genera_to_keep)) %>% #put abundances of genera to keep in a list 
  mutate(community = pmap(., transform_to_taxonomic, cutoff = cutoff_frac)) %>%#When a taxon has less than 5% rel. freq. we put it to zero
  unnest_wider(community) %>%
  filter(community_size>1, period != "1995 - 2000") %>%#delete cases where comm size is only 1 genus because very sad
  select(-n) #ditch raw abundances

# COMPUTE AVAILABLE AND USED RANGE ---------------
available_range <- compute_range_per_combo(data=joined_data_selection, 
                                           combo_vars=c("period", "stress"),
                                           range_vars = chemistry_to_keep,
                                           range_name = "available") %>%
  select(-n)

used_range <- compute_range_per_combo(data=joined_data_selection, 
                                      combo_vars = c("period", "stress", "community"),
                                      range_vars = chemistry_to_keep,
                                      range_name = "used") %>%
  select(-n)
  
used_vs_available_range <- used_range %>%
  left_join(available_range, by = c("period", "stress", "parameter")) %>% #add env variation for corresponding parameter and period
  mutate(parameter = as_factor(gsub("chem_", "", parameter, #nicer names
                                    fixed = TRUE))) %>%
  mutate(parameter = recode_factor(parameter, #nicer names
                                   ammonia_n = "NH4-N", nitrate_n = "NO3-N", ph = "pH")) %>%
  mutate(stress = recode_factor(stress, #nicer names
                                   low = "reference", high = "altered")) %>%
  mutate(community_nr = match(used_range$community, unique(used_range$community))) %>% #give a nr to ID the communities
  mutate(type_of_parameter = ifelse(parameter%in%c("NH4-N", "NO3-N", "silica"), "nutrients", "others")) 
## used vs. available ranges (everything is below id line, av. range almost never entirely used)--------
ggplot(used_vs_available_range) + 
  theme_bw() + 
  theme(panel.grid = element_blank()) +
  scale_colour_manual(values=c(cbPalette)) + 
  geom_boxplot(aes(x=parameter, y=log10(available)-log10(used), col=period)) +
  facet_grid(cols=vars(stress), scales="free") + 
  geom_hline(yintercept = 0, lty="dashed") + 
  labs(x="chemistry variable",
       y=expression(paste(log[10],"(",range[available],"/",range[used],")"))) +
  theme(axis.text.x = element_text(angle = 90))

## used ranges in high vs low stress ------------
used_vs_available_range_main_plot <- used_vs_available_range %>%
  pivot_wider(names_from="stress", values_from = c("used", "available")) %>%
  filter(!is.na(used_reference), !is.na(used_altered), !is.na(available_reference), !is.na(available_altered)) %>%
  mutate(reference_altered = log10(used_reference) - log10(used_altered))

ggplot(used_vs_available_range_main_plot) + 
  theme_bw() + 
  theme(panel.grid = element_blank()) +
  scale_colour_manual(values=cbPalette) + 
  geom_boxplot(aes(x=parameter, y=reference_altered, col=period)) +
  #facet_wrap(vars(period), ncol=4, scales="free") + 
  labs(x="chemistry variable",
     y=expression(paste(log[10],"(",range[reference],"/",range[altered],"), used"))) + 
  theme(axis.text.x = element_text(angle = 90)) + 
  geom_hline(yintercept = 0, lty="dashed")

lm(reference_altered ~ period + parameter , 
            data = used_vs_available_range_main_plot %>%
              mutate(reference_altered = log10(used_reference / used_altered))) %>%
  summary %>% coef 

## check what are these communities (List the genera composition per community)--------
invert_chem_ID <- used_vs_available_range_main_plot %>%
  separate(community, sep="_", into=genera_to_keep, convert=T) %>% #split again genera
  mutate(community_nr = match(used_vs_available_range_main_plot$community_nr, unique(used_vs_available_range_main_plot$community_nr))) %>%
  group_by(community_nr) %>%
  summarise(across(all_of(genera_to_keep), ~ mean(.x))) %>%#compute average membership; just a trick to get rid of multiple lines per community
  select(-3) %>%#kick out the pedigree one
  select(-where(~ is.numeric(.x) && sum(.x)==0)) %>%#get rid of genera w zero density
  pivot_longer(cols=any_of(genera_to_keep)) %>%
  filter(value==1) %>%
  select(-value) 
