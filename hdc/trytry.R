tryCatch({
for (i in 1:3)
{
        cat(i,"\n")
        if (i  == 1) {stop("trying tryCatch")}
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n");next})
