source(file = "./DDelivery_vfinal.R")

raw.data <- import_96wellplates(filepath = "./Raw Data/", 
                                pattern = "Endpoint",
                                remove = 1) 

# ----------------------------------------------------------------
#                  Platemap Annotations
# ----------------------------------------------------------------

polymers <- annotate(wellgrid = c("A1:G1, A5:G5, A9:G9",
                                  "A2:G2, A6:G6, A10:G10", 
                                  "A3:G3, A7:G7, A11:G11", 
                                  "A4:G4, A8:G8, A12:G12"),
                     values = c("dextran", "alpha", "beta", "gamma"),
                     varname = "polymers")
  
crosslinker <- annotate(wellgrid = c("A1:C12", "E1:G12"),
                        values = c("prepolymer", "crosslinked"),
                        varname = "crosslinker")
  
wellcontent <- annotate(wellgrid = c("A1:C12, E1:G12", "H1:H11"),
                        values = c("experimental", "controls"),
                        varname = "wellcontent")

dayplate <- annotate(wellgrid = c("A1:H4", "A5:H8", "A9:H12"),
                     values = c(1,2,3),
                     varname = "dayplate")


annotated_df <- multi_join(polymers, 
                           crosslinker, 
                           wellcontent, 
                           dayplate, 
                           raw.data,  
                           by = "well") %>%
  get_time()

temperature <- annotate(wellgrid = c("B1:C12, F1:G12, H2:H3, H6:H7, H10:H11",
                                     "A1:A12, E1:E12, H1:H1, H5:H5, H9:H9"),
                        values = c(25, 37),
                        varname = "Temperature")

shaking <- annotate(wellgrid = c("C1:C12, G1:G12, H4:H4, H8:H8, H11:H11",
                                 "A1:B12, E1:F12, H1:H2, H5:H6, H9:H10"),
                    values = c(0, 100),
                    varname = "shaking")

annotated_df1 <- annotated_df %>%
  filter(time < 9) %>%
  mutate(shaking = 0) %>%
  mutate(Temperature = 25) 

annotated_df2 <- annotated_df %>%
  filter(time >= 9) %>%
  multi_join(.,temperature, shaking, by = "well") %>%
  bind_rows(.,annotated_df1) %>%
  


experimental.df <- filter(annotated_df, wellcontent == "experimental")

control.df <- filter(annotated_df, wellcontent == "controls")

# ------------------------------------------------------------------------------
#                      Calculations
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
#                     Plotting
# ------------------------------------------------------------------------------

library(ggplot2)

annotated_df2 %>% 
  filter(!is.na(Temperature)) %>%
  filter(wellcontent != "controls") %>% 
  ggplot(., aes(x = time, y = absorbance))+
  geom_point(aes(color = polymers))+
  facet_wrap(shaking~Temperature*crosslinker, nrow = 2)+
  xlab("Time (Days)")+
  theme_bw()

annotated_df2 %>%
  filter(!is.na(absorbance), !is.na(crosslinker)) %>% 
  group_by(Temperature, shaking, time) %>% 
  arrange(time) %>% 
  ggplot(., aes(x = absorbance))+
  geom_histogram(aes(fill = crosslinker), alpha =0.3)+
  facet_wrap(Temperature~shaking)+
  theme_bw()

