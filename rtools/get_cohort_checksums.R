library(dplyr, warn.conflicts=F)
library(tidyr)
library(XML)
library(stringr)
library(purrr)
pattern = "End of step ([0-9]+) of task id ([0-9]+). Checksum = ([0-9.]+)"
parfile = "initparfile.xml"

xml.pars <- xmlToList(parfile)
sp <- xml.pars$sp_name
nb_ages <- sum(as.numeric(unlist(strsplit(xml.pars$nb_cohort_life_stage[[sp]], split=" "))))

get_tstart <- function(task_id){
  if (task_id >= nb_ages){
    tstart_cohort = task_id-nb_ages+1;
  }else{
    tstart_cohort = 0;
  }
  return(tstart_cohort)
}
get_agestart <- function(task_id){
  if (task_id >= nb_ages){
    age_start = 0;
  }else{
    age_start = nb_ages - task_id - 1;
  }
  return(age_start)
}

data <- do.call(rbind, lapply(dir("./", pattern="log_taskfunc"), function(file){
  out <- readLines(file)
  if (length(out)>0){
    out <- out[grepl("Checksum", out)]
    do.call(rbind, lapply(1:length(out), function(i){
      str_match(out[i], pattern)[2:4] %>% 
        t %>% 
        as.data.frame %>% 
        mutate(across(c(V1, V2), ~as.numeric(.x)))
    }))
  }
})) %>% 
  set_names(c("step", "task_id", "checksum")) %>% 
  mutate(t_count = pmap_int(list(step, task_id), .f=function(step, task_id) return(step - get_agestart(task_id) + get_tstart(task_id))),
         age=step) %>% 
  arrange(t_count, age) %>% 
  dplyr::select(t_count, age, checksum) %>% 
  pivot_wider(id_cols="t_count", names_from = "age", values_from="checksum") %>% 
  dplyr::select(-t_count) %>% 
  as.data.frame
print(data)
cat("Writing output to file 'checksum'\n")
data  %>% 
  write.table(file="checksums", row.names=F, quote=F, col.names=F)
