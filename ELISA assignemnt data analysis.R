
## ELISA assignment data analysis 
## Date: 15 Decemebr 2023
## Author: Berhe Etsay 
## Contact: berhe.etsay@gmail.com 



# load packages --------------------------------------------------------------

pacman::p_load(rio, 
               here, 
               tidyverse)


# import data with the raw reading---------------------------------------------

patient_1_raw <- import(("Appendix 2.xlsx"), sheet = "Patient 1") %>% 
  mutate(patient = "Patient 1") %>% 
  mutate(dilution = c(0.10, 0.02, 0.004, 0.0008, 0.00016,0.000032,0.0000064,0))
patient_2_raw <- import(("Appendix 2.xlsx"), sheet = "Patient 2") %>%
  mutate(patient = "Patient 2") %>% 
  mutate(dilution = c(0.10, 0.02, 0.004, 0.0008, 0.00016,0.000032,0.0000064,0))
patient_3_raw <- import(("Appendix 2.xlsx"), sheet = "Patient 3") %>%
  mutate(patient = "Patient 3") %>% 
  mutate(dilution = c(0.10, 0.02, 0.004, 0.0008, 0.00016,0.000032,0.0000064,0))
patient_4_raw <- import(("Appendix 2.xlsx"), sheet = "Patient 4") %>%
  mutate(patient = "Patient 4") %>% 
  mutate(dilution = c(0.10, 0.02, 0.004, 0.0008, 0.00016,0.000032,0.0000064,0))
patient_5_raw <- import(("Appendix 2.xlsx"), sheet = "Patient 5") %>%
  mutate(patient = "Patient 5") %>% 
  mutate(dilution = c(0.10, 0.02, 0.004, 0.0008, 0.00016,0.000032,0.0000064,0))

# Binding raws----------------------------------------------------------------- 

elisa_raw <- bind_rows(patient_1_raw, 
                       patient_2_raw, 
                       patient_3_raw, 
                       patient_4_raw, 
                       patient_5_raw) %>% 
  rename(Ab1_R1 = 3,
         Ab1_R2 = 4,
         Ab2_R1 = 5,
         Ab2_R2 = 6,
         Ab3_R1 = 7,
         Ab3_R2 = 8,
         Ab4_R1 = 9,
         Ab4_R2 = 10,
        `PBS (-ve control)` = 2,
         plate = 1) %>% 
  mutate(dilution = format(dilution,scientific = FALSE)) %>% 
  mutate(dilution = as.numeric(dilution)) %>% 
  rowwise() %>% 
  mutate(
    Average_Ab1 = mean(Ab1_R1,Ab1_R2, na.rm = T),
    Average_Ab2 = mean(Ab2_R1,Ab2_R2, na.rm = T),
    Average_Ab3 = mean(Ab3_R1,Ab3_R2, na.rm = T),
    Average_Ab4 = mean(Ab4_R1,Ab4_R2, na.rm = T)) %>% 
  select(Average_Ab1,
         Average_Ab2,
         Average_Ab3,
         Average_Ab4,
         `PBS (-ve control)`,
         dilution,
         patient 
        )


# Create the long data set  -----------------------------------------------


elisa <- elisa_raw %>% 
  mutate(dilution = log(dilution)) %>% 
  pivot_longer(
    cols = c(Average_Ab1:`PBS (-ve control)`),
    names_to = "Antibody",
    values_to = "Absorbtion"
  ) %>% 
  mutate(Absorbtion = log(Absorbtion)) %>% 
  mutate(Antibody = recode(Antibody,
                           "Average_Ab1" = "Antibody 1",
                           "Average_Ab2" = "Antibody 2", 
                           "Average_Ab3" = "Antibody 3", 
                           "Average_Ab4" = "Antibody 4" 
  )) %>% 
  mutate(identity_combined = paste(patient, Antibody))

# create the background noise for the standard (row-H) 

Background_noise <- elisa %>% 
  filter(dilution == -Inf & Antibody != "PBS (-ve control)") %>% 
  mutate(identity_combined = paste(patient, Antibody))



# plot the concentration------------------------------------------------------
elisa %>% 
  filter(dilution != -Inf) %>% 
ggplot(aes(x = dilution,
           y = Absorbtion,
           color = Antibody)) +
  geom_point() +
  geom_line(size = 1.0) +
  scale_x_continuous(breaks = c(-11.9, -10.3, -8.7, -7.1, -5.5,-3.9,-2.3)) +
  facet_wrap(~patient,
             scales = "free",
             ncol = 2) +
  labs(y = "The log of average spectrophotmeter reading for each antibody",
       x = "The log of dilution (from dilution factor 1/156,250 to 1/10)",
       title = "Indirect sandwich ELISA test of four different antibodies to diagnose DENV-3 from patient sera sample",
       subtitle = "**Scales are in logarithmic scale",
       fill = "Antibody",
       caption = "The photospectometer reading and the dilution factor are changed in to base log values 
       from -11.9 (1/156250) to -2.3 (1/10) in the x axis, and from -1.6 (0.1946) to 0.5 (1.5995) in the y axis") +
  theme(axis.title = element_text(size = 12, face = "italic", color = "blue"),
        plot.title = element_text(size = 12, face = "bold", color = "red"),
        plot.subtitle = element_text(size = 12, face = 'italic', color = "blue"),
        strip.placement = "outside", 
        strip.text=element_text(face ="italic", color = "red", size = 10),
        plot.caption = element_text(face = "italic", color = "red"))


# Display background noise ------------------------------------------------

elisa %>% 
  filter(dilution != -Inf) %>% 
  ggplot(aes(x = dilution,
             y = Absorbtion,
             color = Antibody)) +
  geom_point() +
  geom_line(size = 1.0) +
  scale_x_continuous(breaks = c(-11.9, -10.3, -8.7, -7.1, -5.5,-3.9,-2.3)) +
  facet_wrap(~patient,
             scales = "free",
             ncol = 2) +
  labs(y = "The log of average spectrophotmeter reading for each antibody",
       x = "The log of dilution (from dilution factor 1/156,250 to 1/10)",
       title = "Indirect sandwich ELISA test of four different antibodies to diagnose DENV-3 from patient sera sample",
       subtitle = "**Scales are in natural logarithmic scale **not base10**",
       fill = "Antibody",
       caption = "The photospectometer reading and the dilution factor are changed in to base log values 
       from -11.9 (1/156250) to -2.3 (1/10) in the x axis, and from -1.6 (0.1946) to 0.5 (1.5995) in the y axis") +
  theme(axis.title = element_text(size = 12, face = "italic", color = "blue"),
        plot.title = element_text(size = 12, face = "bold", color = "red"),
        plot.subtitle = element_text(size = 12, face = 'italic', color = "blue"),
        strip.placement = "outside", 
        strip.text = element_text(face = "italic", color = "red", size = 10),
        plot.caption = element_text(face = "italic", color = "red"))


# Display background noise ------------------------------------------------

elisa %>% 
  filter(dilution != -Inf & Antibody != "PBS (-ve control)") %>% 
  ggplot(aes(x = dilution,
             y = Absorbtion,
             color = Antibody)) +
  geom_point() +
  geom_line(size = 1.0) +
  geom_hline(data = Background_noise,
            aes(yintercept = Absorbtion,lty = "Averge reading for \n *PBS + Capture Ab + Detection Ab*")) +
  geom_hline(data = patients_naive_sera,
             aes(yintercept = mean_ab_naive_s, lty = "Average reading for \n*Capture Ab + Naiva sera*" )) +
  scale_x_continuous(breaks = c(-11.9, -10.3, -8.7, -7.1, -5.5,-3.9,-2.3)) +
  facet_grid(Antibody~patient) +
  labs(y = "The log of average spectrophotmeter reading for each antibody",
       x = "The log of dilution (from dilution factor 1/156,250 to 1/10)",
       title = "Indirect sandwich ELISA test of four different antibodies to diagnose DENV-3 from patient sera sample",
       subtitle = "**Scales are in natural logarithmic scale **not base10**",
       fill = "Antibody",
       caption = "** The horizontal dashed lines are the balnk wells reading
       where there is no antigen added but there is PBS + capture antibody and
       secondary detection antibodies") +
  theme(axis.title = element_text(size = 12, face = "italic", color = "blue"),
        plot.title = element_text(size = 12, face = "bold", color = "red"),
        plot.subtitle = element_text(size = 12, face = 'italic', color = "blue"),
        strip.placement = "outside", 
        strip.text=element_text(face ="italic", color = "red", size = 10),
        plot.caption = element_text(face = "italic", color = "red"))



