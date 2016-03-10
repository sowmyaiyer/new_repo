values <- as.numeric(scan(file="stdin", quiet=TRUE))
q <- quantile(values,prob=seq(0,1,0.25))
cat(q)
