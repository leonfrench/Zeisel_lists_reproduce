library(xlsx)
library(homologene)
library(readxl)
library(reshape2)
library(readr)
library(here)
library(dplyr)
library(tidyr)
library(magrittr)
library(conflicted)
conflict_prefer("first", "dplyr")
conflict_prefer("count", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("flatten", "purrr")

#read in Ziesel supplement file
mouse_supplement <- read_xlsx(here("data","Zeisel 2015 mouse original", "aaa1934_tables1.xlsx"))
mouse_supplement %<>% mutate(row_num = as.character(1:n()))
mouse_supplement <- melt(mouse_supplement, id.vars = "row_num") %>% tibble()
mouse_supplement %<>% select(-row_num)
mouse_supplement %<>% filter(value != "All genes")
mouse_supplement %<>% rename(cell_type = variable, gene_symbol = value)
mouse_summary <- mouse_supplement %>% group_by(cell_type) %>% count() %>% rename(mouse_n = n)
mouse_summary


#convert to human, this is the only part that is dynamic, everything else is static. However, gene symbols will be out of date.
human_mapping <- mouse2human(mouse_supplement %>% pull(gene_symbol)) %>% tibble()
human_converted <- inner_join(mouse_supplement, human_mapping %>% rename(gene_symbol = mouseGene))
human_converted %<>% select(gene_symbol = humanGene, cell_type) %>% distinct()

#old_interneuron <- read_csv("/Users/leon/Downloads/UpdatedCellTypeLists/old Zeisel lists/Zeisel.interneuron.txt", col_names = F)
#old_interneuron %<>% rename(gene_symbol = X1) %>% mutate(cell_type = "Old interneuron list")
#human_converted <- bind_rows(human_converted, old_interneuron)

human_converted_summary <- human_converted %>% group_by(cell_type) %>% count() %>% rename(human_n = n)
human_converted_summary

#filter for Allen 6 - stage 1, this file is from figshare
allen_6 <- read_tsv(here("data","Allen 6 correlation Frontiers","AllenHBA_DK_ExpressionMatrix.tsv"))
allen_6 %<>% rename(gene_symbol = ...1)



#take a look at the interneurons
intersect(allen_6$gene_symbol, human_converted %>% filter(cell_type == "Interneuron") %>% pull(gene_symbol)) %>% length()
human_converted %>% filter(cell_type == "Interneuron") %>% pull(gene_symbol) %>% length()
allen_6 %<>% filter(`Average donor correlation to median` > 0.446)


human_stage1 <- human_converted %>% filter(gene_symbol %in% allen_6$gene_symbol)
human_stage1_summary <- human_stage1 %>% group_by(cell_type) %>% count() %>% rename(stage1_n = n)
human_stage1_summary

#filter for BrainSpan correlations - stage 2
brainspan_cor <- read_csv(here("data","BrainSpan correlations","BrainSpanCorrelations.csv"))
brainspan_cor %<>% filter(pvalue < 0.05)
human_stage2 <- human_stage1 %>% filter(gene_symbol %in% brainspan_cor$GeneSymbol)
human_stage2_summary <- human_stage2 %>% group_by(cell_type) %>% count() %>% rename(stage2_n = n)
human_stage2_summary

#test against supplement from Shin et al. 2015 cortical thickness paper
#these are the gene sets after two stage filtering steps
reference_from_paper <- read_tsv(here("data","From Figshare of thickness paper","Reference_Consistent_Genes_ObtainedBy2StageFiltering.tsv"))
reference_from_paper_summary <- reference_from_paper %>% group_by(CellType) %>% count() %>% filter(!is.na(CellType)) %>% rename(cell_type = CellType, n_from_Shin_2015 = n)
reference_from_paper_summary %<>% mutate(cell_type = gsub("[.]", " ", cell_type))
inner_join(human_stage2_summary, reference_from_paper_summary) #matches

#write out full table
joined_summaries <- inner_join(inner_join(inner_join(inner_join(mouse_summary, human_converted_summary) ,human_stage1_summary), human_stage2_summary), reference_from_paper_summary)
joined_summaries

dir.create(here("results"), showWarnings = F)
joined_summaries %>% write_csv(here("results", "combined_counts.csv"))

#write out lists
mouse_supplement %>% write_csv(here("results", "original_mouse_lists.csv"))
human_converted %>% write_csv(here("results", "human_converted_lists.csv"))
human_stage1 %>% write_csv(here("results", "Stage1_lists.csv"))
human_stage2 %>% write_csv(here("results", "Stage2_lists.csv"))

write.xlsx(mouse_supplement, here("results", "original_mouse_lists.xlsx"), sheetName = "Sheet1", col.names = TRUE, row.names = T, append = FALSE)
write.xlsx(human_converted, here("results", "human_converted_lists.xlsx"), sheetName = "Sheet1", col.names = TRUE, row.names = T, append = FALSE)
write.xlsx(human_stage1, here("results", "Stage1_lists.xlsx"), sheetName = "Sheet1", col.names = TRUE, row.names = T, append = F)
write.xlsx(human_stage2, here("results", "Stage2_lists.xlsx"), sheetName = "Sheet1", col.names = TRUE, row.names = T, append = F)
