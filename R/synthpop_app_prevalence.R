library(tidyverse)

app_matrix <- read.csv("~/Downloads/empop_ego_app_matrix_binary.csv")
demos <- read.csv("~/Downloads/empop_egoid_demo_def.csv")

app_long <- app_matrix %>%
  pivot_longer(starts_with("a"),
               names_to = "app",
               values_to = "use") %>%
  filter(use == 1) %>%
  rename(egoid = X) %>%
  left_join(demos, by = "egoid") %>%
  select(app, demo8, use) %>%
  group_by(app, demo8) %>%
  summarize(num_using = sum(use))

demo_sizes <- demos %>%
  group_by(demo8) %>%
  summarize(demo_total = n())

demo_labels <- demos %>%
  group_by(demo8) %>%
  slice(1) %>%
  mutate(stratum = paste(race_ethnicity, age, sep = "_")) %>%
  select(demo8, stratum)

app_merge <- app_long %>%
  left_join(demo_sizes, by = "demo8") %>%
  mutate(prev = num_using/demo_total) %>%
  left_join(demo_labels, by = "demo8") %>%
  select(app, prev, stratum) %>%
  pivot_wider(id_cols = app,
              names_from = stratum,
              values_from = prev,
              values_fill = 0)

write.csv(app_merge, "~/Desktop/chiSTIG_HPC/R/synthpop_app_prevalence.csv")
