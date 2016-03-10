q <- p.adjust(as.numeric(scan(commandArgs(TRUE))), method="BH")
write(q,
