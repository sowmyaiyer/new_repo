BEGIN{ FS="\t"; OFS="\t"}
{
	if ($7 == "-") { tss = $4} else {tss = $5}
	split($9,arr,";")
	for (i = 1; i <= length(arr); i ++)
	{
		if (index(arr[i],"FPKM ") != 0)
		{	
			split(arr[i],fpkm," ")
			gsub("\"","",fpkm[2])
			print ($1,tss,tss,fpkm[2])
		}
	}
}
