# library(tximeta)
library(tximport)
library(dplyr)
library(reshape2)

# files = file.path("SRR04_output/alevin/quants_mat.gz")
# file.exists(files)
# coldata = data.frame(files, names="SRR04_matrix", stringsAsFactors=FALSE)
# file04_matrix = tximeta(coldata)
srr04_file = file.path("SRR04_output/alevin/quants_mat.gz")
srr05_file = file.path("SRR05_output/alevin/quants_mat.gz")
srr05_file = file.path("SRR05_output/alevin/quants_mat.gz")

srr04 = tximport(files=srr04_file, type="alevin")
srr05 = tximport(files=srr05_file, type="alevin")
srr06 = tximport(files=srr06_file, type="alevin")

srr04_matrix_only = srr04$counts
srr05_matrix_only = srr05$counts
srr06_matrix_only = srr06$counts


writeMM(srr04_matrix_only, "srr04.mtx", col.names = TRUE, row.names = TRUE)
writeMM(srr04_matrix_only, "srr04.mtx", col.names = TRUE, row.names = TRUE)
writeMM(srr04_matrix_only, "srr04.mtx", col.names = TRUE, row.names = TRUE)

readssr04mtx = readMM("srr04.mtx")