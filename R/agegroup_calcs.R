age_race %>% pivot_wider(names_from = "Var2", values_from = "Freq",
                         names_prefix = "agegrp") %>%
  mutate(exp1 = ((agegrp1 + agegrp2)/14) * 5,
         exp2 = ((agegrp1 + agegrp2)/14) * 9,
         pct1 = agegrp1/(agegrp1 + agegrp2),
         pct2 = agegrp2/(agegrp1 + agegrp2),
         exppct1 = exp1/(agegrp1 + agegrp2),
         exppct2 = exp2/(agegrp1 + agegrp2)
  )


netstats <- readRDS("~/Desktop/chiSTIG_hpc/data/intermediate/estimates/netstats-local.rds")


netstats_table <- as.data.frame(table(netstats$attr$race, netstats$attr$age.grp)) %>% pivot_wider(names_from = "Var2", values_from = "Freq",
                                                                                                 names_prefix = "agegrp") %>%
  mutate(exp1 = ((agegrp1 + agegrp2)/14) * 5,
         exp2 = ((agegrp1 + agegrp2)/14) * 9,
         pct1 = agegrp1/(agegrp1 + agegrp2),
         pct2 = agegrp2/(agegrp1 + agegrp2),
         exppct1 = exp1/(agegrp1 + agegrp2),
         exppct2 = exp2/(agegrp1 + agegrp2)
  ) %>% rename("race" = Var1)

netstats_table

demo_counts %>% group_by(race_ethnicity) %>%
  summarize(per_year = sum(population)/14)
