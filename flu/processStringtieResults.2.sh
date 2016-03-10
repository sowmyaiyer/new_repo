for time in {"6h","12h","24h","48h","mock12h"}
do
	echo $time
	sed 's/\"//g' /home/si14w/gnearline/flu/stringtie_out/stringtie_out_${time}.gtf | sed 's/;$//g' | awk -F"\t" '
BEGIN {OFS="\t"; print "gene","transcript","length","TPM"}
{
	if ($3 == "transcript") {
	split($NF,arr,"; ");
	len = $5-$4+1
	split(arr[1], strg_gene_nm_arr, " ");
	split(arr[2], strg_transcript_nm_arr, " ");
	stringtie_gene_id = strg_gene_nm_arr[2];
	stringtie_transcript_id = strg_transcript_nm_arr[2];
	split(arr[4], ref_gene_nm_arr, " ");
	split(arr[3], ref_transcript_nm_arr," ");
	split(arr[length(arr)], tpm_arr," ");
	tpm = tpm_arr[2];
	if (ref_gene_nm_arr[1] == "ref_gene_id")
		geneid = ref_gene_nm_arr[2];
	else
		geneid = stringtie_gene_id;
	if (ref_transcript_nm_arr[1] == "reference_id")
		transcriptid = ref_transcript_nm_arr[2];
	else 
		transcriptid = stringtie_transcript_id;
	print geneid,transcriptid,len,tpm
	}
}'  > ../txt/tpm_${time}.txt
done
