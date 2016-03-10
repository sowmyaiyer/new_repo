for (iter in 1:500000)
{
	letter <- character(10)
	len <- sample(c(9,10,11,12,13,14,15), prob=c(0.1,0.2,0.3,0.2,0.1,0.05,0.025), size=1)
	for (i in 1:len)
	{
		letter[i] <- sample(c("A","T","C","G"), prob=c(0.25,0.25,0.25,0.25), size=1)
	}
	cat (paste(">random10mer_",iter,sep="","\n"),paste(letter, collapse=""),sep="","\n")
}
