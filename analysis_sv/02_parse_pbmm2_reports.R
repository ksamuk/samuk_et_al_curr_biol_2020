# parse pbmm2 output

library("tidyverse")

pb_files <- list.files("analysis_sv/pacbio_out", full.names = TRUE)

pb_file <- pb_files[1]


parse_pbmm2_reports <- function(pb_file){
  
  message(pb_file)
  
  pb_lines <- read_lines(pb_file)
  
  start_line <- grep("Mapped Reads:", pb_lines)
  end_line <- grep("Peak RSS:", pb_lines)
  
  if(length(start_line) == 1 & length(end_line) == 1){
    
    pb_lines <- pb_lines[c(2,start_line:end_line)]
    pb_lines <- gsub(".*-\\|- ", "", pb_lines)
    
    barcode <- gsub("sample ", "", pb_lines[1])
    
    pb_lines <-  pb_lines[-1]
    
    data.frame(barcode, pb_lines) %>%
      separate(pb_lines, into = c("stat", "value"), sep = ": ") 
    
    
  }
  
}

pb_reports <- lapply(pb_files, parse_pbmm2_reports)
pb_reports <- bind_rows(pb_reports)

pb_reports %>%
  filter(stat == "Mapped Reads") %>%
  ggplot(aes(x = barcode, y = as.numeric(value)))+
  geom_point(size = 3)+
  facet_wrap(~stat, scales = "free_y")+
  ylab("Mapped Bases")

pb_reports %>%
  filter(stat == "Mapped Bases") %>% 
  mutate(value = as.numeric(value)) %>%
  mutate(barcode = fct_reorder(barcode, value)) %>%
  ggplot(aes(x = barcode, y = value/1000000000))+
  geom_point(size = 3)+
  facet_wrap(~stat, scales = "free_y")+
  ylab("Mapped Bases (GB)")+
  xlab("Barcode Index") +
  ggtitle("Bases mapped to the D. pseudoobscura reference genome")

