# Download all database in ExmRes.1.ANNOVAR_download_databases.sh
for i in $(seq 1 24)
do
	bash ExmRes.1.ANNOVAR_download_databases.b38.sh $i
done
