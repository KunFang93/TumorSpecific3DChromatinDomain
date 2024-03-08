#$1 is resolution
half_res=$(($1/2))
for i in NT1 NT2 PT1 PT2 PT3 PT4 PT5 RT1 RT2 RT3 RT4 RT5
do
    echo $i
    mkdir -p "$1"/"$i"
    # hicpro2fithic
    python /data/kun/Softwares/fithic/fithic/utils/HiCPro2FitHiC.py -i /data/kun/Lava/hicpro_brca_tissue_final/hicpro_matrix/raw/"$i"_"$1".matrix -b /data/kun/Lava/hicpro_brca_tissue_final/hicpro_matrix/raw/"$i"_"$1"_abs.bed -s /data/kun/Lava/hicpro_brca_tissue_final/hicpro_matrix/iced/"$i"_"$1"_iced.matrix.biases -o "$1"/"$i" -r "$1"
    # run fithic 
    L_para=$(($1*2))
    python /data/kun/Softwares/fithic/fithic/fithic.py -i /data/kun/Lava/SLG-DEGs/fithic/tissue/"$1"/"$i"/fithic.interactionCounts.gz -f /data/kun/Lava/SLG-DEGs/fithic/tissue/"$1"/"$i"/fithic.fragmentMappability.gz -t /data/kun/Lava/SLG-DEGs/fithic/tissue/"$1"/"$i"/fithic.biases.gz -r "$1" -L $L_para -x intraOnly -U 1000000 -v -l "$i"_fithic -o /data/kun/Lava/SLG-DEGs/fithic/tissue/"$1"/"$i"
    # merge 
    sh /data/kun/Softwares/fithic/fithic/utils/merge-filter-parallelized.sh /data/kun/Lava/SLG-DEGs/fithic/tissue/"$1"/"$i"/"$i"_fithic.spline_pass1.res"$1".significances.txt.gz $1 /data/kun/Lava/SLG-DEGs/fithic/tissue/"$1"/"$i"/"$i"_fithic-"$1"-merge 0.01 /data/kun/Softwares/fithic/fithic/utils/
    # post merge
    zcat /data/kun/Lava/SLG-DEGs/fithic/tissue/"$1"/"$i"/"$i"_fithic-"$1"-merge/*/subset_fithic_*.gz |sort -k1,1V -k2,2n -k3,3n - > /data/kun/Lava/SLG-DEGs/fithic/tissue/"$1"/merged_results/"$i"_fithic_merged.txt
    # post merge expand
    awk -v res=$half_res 'BEGIN{OFS="\t"}{print $1,$2-res,$2+res,$3,$4-res,$4+res,$5}' "$1"/merged_results/"$i"_fithic_merged.txt > "$1"/merged_results/"$i"_fithic_merged.expand.bed
done