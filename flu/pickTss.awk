{
	split($5, intervalIds, ",")
	split($6, scores, ",")
	num=length(intervalIds)
	if (num == 1)
	{
		printf("%s\t%d\t%d\t%s\t%.2f\t%s\n",$1,$2,$3,$5,$6,$4)
	}
	else {
		argmax=-1
		max=-1
		for (i = 1 ; i <= num; i ++)
		{
			if (scores[i] > max)
			{
				max=scores[i]
				argmax=i
			}
		}
		split(intervalIds[argmax], interval,"@")
		printf("%s\t%d\t%d\t%s\t%.2f\t%s\n",interval[1],interval[2],interval[3],intervalIds[argmax],max,$4)
	}
}
