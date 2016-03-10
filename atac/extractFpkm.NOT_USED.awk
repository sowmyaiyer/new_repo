BEGIN{ FS="\t"; OFS="\t"}
{
	split($14,arr,";")
	for (i = 1; i <= length(arr); i ++)
	{
		if (index(arr[i],"FPKM ") != 0)
		{	
			split(arr[i],fpkm," ")
			gsub("\"","",fpkm[2])
			print (type,fpkm[2])
		}
	}
}
