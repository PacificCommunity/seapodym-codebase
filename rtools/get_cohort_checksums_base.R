library(XML)

parfile = "initparfile.xml"

xml.pars <- xmlToList(parfile)
sp <- xml.pars$sp_name
nb_ages <- sum(as.numeric(unlist(strsplit(xml.pars$nb_cohort_life_stage[[sp]], split=" "))))
 
files <- list.files(pattern="log_taskfunc", full.names=TRUE)
lines <- unlist(lapply(files, readLines))
clines <- lines[grep("Checksum =", lines)]

pattern <- "End of step ([0-9]+) of task id ([0-9]+)\\. Checksum = ([0-9.]+)"
df <- strcapture(pattern, clines, data.frame(age=numeric(), id=numeric(), csum=numeric()))

is.new <- df$id >= nb_ages
df$tcount  <- df$age - ifelse(is.new, 0, nb_ages-df$id-1) + ifelse(is.new, df$id-nb_ages+1, 0)

mat <- tapply(df$csum, list(df$tcount, df$age), max)
def.options <- options(digits = 16)

colnames(mat) <- paste0("checksum", 1:nb_ages-1)
rownames(mat) <- paste0("time", 1:nrow(mat))
print(mat)

cat("Writing output to file 'checksums'\n")
m <- matrix(sprintf("%.16g", mat), nrow = nrow(mat))
write.table(m, file="checksums", row.names=F, col.names=F, quote=F)

options(def.options)

