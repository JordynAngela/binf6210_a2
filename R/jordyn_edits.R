#*************************************************************************************
#Your name: Jordyn Angela Anene Your Student ID:
#Name of peer who was primary author of the script you edited:  Ali Al-Reda El-Zein
#Your GitHub link: https://github.com/JordynAngela
#GitHub link of peer: https://github.com/AliZein-1
#*************************************************************************************
## Load packages -----

rm(list = ls())

library(tidyverse)


## Import data -----

df_trich_s <- read_csv("data/trichinella_spiralis_te.csv")
genus_types <- read_csv("data/2026-10-19_trich_genotypes.csv")
trich_zb <- read_csv("data/trich_z.csv")

## filtering data for DNA transposons for T. zimbabwensis -----

trich_zb$matching_class |> unique()

non_te_z <- c("Simple_repeat", "Low_complexity", "Satellite", "snRNA", "rRNA", "Unknown")

dna_te_z <- c("DNA/PIF-Harbinger", "DNA/TcMar-Tc1", "DNA/MULE-MuDR", 
              "DNA/TcMar-Tc4", "DNA/TcMar-Pogo", "DNA/CMC-Chapaev", 
              "DNA/TcMar-Mariner", "DNA/PiggyBac","RC/Helitron")

trich_z <- trich_zb |> 
  filter(matching_class %notin% non_te_z) |>
  filter(matching_class %in% dna_te_z) 

trich_z$species = c("Trichinella zimbabwensis")

## filtering data for DNA transposons for T. spiralis -----

df_trich_s$matching_class |> unique()

non_te_s <-  c("simple_repeat", "Simple_repeat", "Low_complexity", "Satellite", "rRNA", "snRNA", "ARTEFACT", "Unknown")

dna_te_s <- c("DNA/TcMar-Tc4", "DNA/TcMar-Tc1", "DNA/hAT", 
              "DNA/PIF-Harbinger", "DNA/TcMar-Mariner", "RC/Helitron", "DNA/hAT-Ac", "DNA/MULE-MuDR", "DNA/TcMar-Pogo", "DNA/CMC-Chapaev")

trich_s <- df_trich_s |> 
  filter(matching_class %notin% non_te_s) |> 
  filter(matching_class %in% dna_te_s) 

trich_s$species = c("Trichinella spiralis")


## identifying and filtering the common DNA transposons between the two species -----

trich_s$matching_class |> unique()

trich_z$matching_class |> unique()

dte_common <- c("DNA/PIF-Harbinger", "DNA/TcMar-Tc1", "DNA/MULE-MuDR", 
                "DNA/TcMar-Tc4", "DNA/TcMar-Pogo", "DNA/CMC-Chapaev", 
                "DNA/TcMar-Mariner", "RC/Helitron")
#JA Edit:
#JA# Replacing this to keep all common classes# ts_common <- trich_s |> 
#JA# Replacing this to keep all common classes# group_by(matching_class) |> 
#JA# Replacing this to keep all common classes# select(matching_class, qry_width, perc_insert, species) |> 
#JA# Replacing this to keep all common classes# filter(matching_class %in% dte_common) |> 
#JA# remove the following 2 lines because this actually changes the plot shown 
#JA# (remove this) Don't deduplicate by the measurements# distinct(perc_insert, .keep_all = TRUE) |> 
#JA# (remove this)  Don't deduplicate by the measurements# distinct(qry_width, .keep_all = TRUE)
#JA#tz_common <- trich_z |> 
#JA# (removed) group_by(matching_class) |> 
#JA# (removed) select(matching_class, qry_width, perc_insert, species) |> 
#JA# (removed) filter(matching_class %in% dte_common) |> 
#JA# (removed) distinct(perc_insert, .keep_all = TRUE) |> 
#JA# (removed) distinct(qry_width, .keep_all = TRUE)


#JA# Keep every hit in the common DNA transposon classes
ts_common <- trich_s |>
  filter(matching_class %in% dte_common) |>
  select(
    matching_class, qry_id, qry_start, qry_end, repeat_id,
    in_repeat_start, in_repeat_end, qry_width, perc_insert, species
  )

tz_common <- trich_z |>
  filter(matching_class %in% dte_common) |>
  select(
    matching_class, qry_id, qry_start, qry_end, repeat_id,
    in_repeat_start, in_repeat_end, qry_width, perc_insert, species
  )
#JA# this code keeps all hits for the common classes
#JA# then selects the largest hit or hits per class from the full data.
ts_largest <- ts_common |>
  group_by(matching_class) |>
  slice_max(order_by = qry_width, n = 1, with_ties = TRUE) |>
  ungroup()

tz_largest <- tz_common |>
  group_by(matching_class) |>
  slice_max(order_by = qry_width, n = 1, with_ties = TRUE) |>
  ungroup()
#Ja# the with_ties keep any max width results that tied
## filtering the 5 largest DNA transposons in the two species -----

ts_5 <- ts_common|> 
  top_n(5, qry_width)

tz_5 <- tz_common |> 
  top_n(5, qry_width)

combined_5 <- rbind(ts_5, tz_5)


## filtering the largest transposons for each species -----

ts_largest <- ts_common |> 
  slice(which.max(qry_width))  

tz_largest <- tz_common |> 
  slice(which.max(qry_width))

## creating a combined dataframe for the largest transposons in both species -----

combined_largest <- rbind(tz_largest, ts_largest)

## creating T. zimbabwensis plot -----

p_tz <- ggplot(
  tz_largest,
  aes(x = qry_width, y=perc_insert, colour = matching_class)
) +
  geom_point(size = 5) +
  labs(title = "How Sequence Length Affects Percent Insertion in T. zimbabwensis for Each DNA Transposon Type.",
       subtitle = "Using the Largest Sequence per Transposon Type (In Base Pairs)",
       x = "Sequence Length (bp)",
       y = "Percent Insertion",
       colour = "DNA Transposon Type")

print(p_tz)

ggsave("tz.png", path = "../output")

## creating T. spiralis plot -----
p_ts <- ggplot(
  ts_largest,
  aes(x = qry_width, y=perc_insert, colour = matching_class)
) +
  geom_point(size = 5)+
  labs(title = "How Sequence Length Affects Percent Insertion in T. spiralis for Each DNA Transposon Type.",
       subtitle = "Using the Largest Sequence per Transposon Type (In Base Pairs)",
       x = "Sequence Length (bp)",
       y = "Percent Insertion",
       colour = "DNA Transposon Type")

print(p_ts)   

ggsave("ts.png", path = "../output")

## creating combined plot for top 5 largest sequences per DNA transposon -----

combined_top5 <- ggplot(
  combined_5,
  aes(x =qry_width, y =perc_insert, 
      colour = matching_class, shape = combined_5$species)) +
  geom_point(size = 3) + 
  geom_line(size = 0.5) +
  xlim(0,750) +
  ylim(0,15) +
  labs(title = "The Effects of Sequence Length on Percent Insertion in T. spiralis & T.zimbabwensis for Each DNA Transposon Type.",
       subtitle = "Using the Top 5 Largest DNA Transposon Sequences in Base Pairs.",
       x = "Sequence Length (bp)",
       y = "Percent Insertion",
       colour = "DNA Transposon Type",
       shape = "Species")


print(combined_top5)

ggsave("combined_top5.png", path = "../output")


## creating combined plot for largest transposon length per DNA transposon-----

combined_largest <- ggplot() +
  geom_point(data = combined_largest, 
             aes(x =qry_width, y =perc_insert, 
                 colour = matching_class, shape = combined_largest$species,
                 size = 5)) +
  guides(size = "none") +
  labs(title = "The Effects of Sequence Length on Percent Insertion in T. spiralis & T.zimbabwensis for Each DNA Transposon Type.",
       subtitle = "Using the Largest DNA Transposon Sequences in Base Pairs.",
       x = "Sequence Length (bp)",
       y = "Percent Insertion",
       colour = "DNA Transposon Type",
       shape = "Species")

print(combined_largest)

ggsave("combined_largest.png", path = "../output")

