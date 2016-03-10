for time in {"A","B","C","D","E","F","G","H"}
do
        echo $time
        awk '{ if ($5 > 5) print }' /home/si14w/gnearline/flu/txt/${time}_raw.bed >  /home/si14w/gnearline/flu/txt/${time}_raw.tags_above_5.bed
done
