for samplename in  "NT5" "PT4" "PT5" "RT1_022421" "RT2_022421" 
do
    for binsize in 20000 40000 100000 150000
    do
        for winsize in 5 8 10 13 15
        do
            mkdir -p /data/kun/Lava/TADs_calling/"$samplename"
            echo "$samplename $binsize $winsize" 
            Rscript /data/kun/Lava/process_codes/topdom_rscript.R /data/kun/Lava/hicpro_result_Jin_V_022421/hic_results/matrix/"$samplename"/iced/"$binsize" /data/kun/Lava/TADs_calling/"$samplename"/"$binsize"/"win$winsize" $winsize
        done
    done
done
