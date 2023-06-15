# READ THE COUPLED CHEMISTRY-INVERT DATA ------
joined_data_selection <- readRDS("joined_data_selection.rds") %>%#  joined_data_selection %>%
  select(-n) #ditch raw abundances

# COMPUTE AVAILABLE AND USED RANGE ---------------
chemistry_to_keep <- c("chem_ph", "chem_cond", "chem_turb", "chem_doc",
                       "chem_calcium", "chem_ammonia_n", "chem_nitrate_n", "chem_silica")
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
## used vs. available ranges (indicating that available range is almost never entirely used)--------
ggplot(used_vs_available_range) + 
  theme_bw() + 
  theme(panel.grid = element_blank()) +
  scale_colour_manual(values=c(cbPalette)) + 
  aes(x=log10(available), y=log10(used), 
      col=parameter, pch=period)+
  geom_point() +
  facet_grid(cols=vars(stress), scales="free") + 
  geom_abline(intercept = 0, slope = 1, lty="dashed") + 
  labs(x=expression(paste(log[10],"(",range[available],")")),
       y=expression(paste(log[10],"(",range[used],")"))) 

## used ranges in high vs low stress ------------
used_vs_available_range_main_plot <- used_vs_available_range %>%
  pivot_wider(names_from="stress", values_from = c("used", "available")) %>%
  filter(!is.na(used_reference), !is.na(used_altered), !is.na(available_reference), !is.na(available_altered)) %>%
  mutate(reference_altered = log10(used_reference) - log10(used_altered))

ggplot(used_vs_available_range_main_plot) + 
  aes(x=log10(used_reference), y=log10(used_altered), 
      col=parameter, pch=period)+
  theme_bw() + 
  theme(panel.grid = element_blank()) +
  scale_colour_manual(values=cbPalette) + 
  geom_point() +
  facet_wrap(vars(type_of_parameter), ncol=2, scales="free") + 
  theme(strip.text.x = element_text(size = 15)) +
  labs(x=expression(paste(log[10],"(",range[reference],")")),
       y=expression(paste(log[10],"(",range[altered],")"))) + 
  geom_abline(intercept=0, slope=1, lty="dashed") +
  guides(col = guide_legend(ncol = 2))

lm(reference_altered ~ period + parameter , 
            data = used_vs_available_range_main_plot %>%
              mutate(reference_altered = log10(used_reference / used_altered))) %>%
  summary %>% coef 
