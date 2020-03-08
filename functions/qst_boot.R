library("tidyverse")
library("lme4")

co_counts <- read.table("data/co_counts_manual.txt", h = T)

###########################################
# QST
###########################################

# QST CALCULATION

compute_qst_boot <- function(filtered_df, reps = 100, out_folder = "data/qst_boot"){
  
  dir.create(out_folder)
  
  
  qst_boot <- rep(NA, reps)
  
  pop1 <- filtered_df %>%
    mutate(pop = as.factor(pop), line = as.factor(line)) %>%
    filter(pop == "AFC")
  
  pop2 <- filtered_df %>%
    mutate(pop = as.factor(pop), line = as.factor(line)) %>%
    filter(pop == "MC")
  
  for (i in 1:reps){

    pop1_boot <- pop1  %>%
      sample_frac(1, replace = TRUE) 
    
    pop2_boot <- pop2  %>%
      sample_frac(1, replace = TRUE)
    
    filtered_df_boot <- bind_rows(pop1_boot, pop2_boot)
    
    #var_glm <- blmer(crossovers ~ (1|pop/line), data = filtered_df_boot,  control = lmerControl(optimizer ="Nelder_Mead"))
    var_glm <- lmer(crossovers ~ (1|pop/line), data = filtered_df_boot,  control = lmerControl(optimizer ="Nelder_Mead"))
    #var_glm <- glmer(crossovers ~ (1|pop/line), family = "poisson", data = filtered_df)
    
    if(!isSingular(var_glm)){
      
      var_comp <- VarCorr(var_glm) %>% data.frame
      var_pop <- var_comp$sdcor[2]
      var_cross <- var_comp$sdcor[1]
      
      qst_boot[i] <- var_pop^2 / (var_pop^2 + var_cross^2)
      
    }
    
    
  }
  
  file_uuid <- stringi::stri_rand_strings(1, 40)
  file_name <- paste0(out_folder, "/qst_boot_rep_", file_uuid, ".txt")
  qst_boot_df <- data.frame(qst_boot = qst_boot)
  
  write.table(qst_boot_df, file_name, col.names = TRUE, row.names = FALSE, quote = FALSE)
  
}

for (i in 1:100){
  
  compute_qst_boot(co_counts, reps = 100)
  
}



