python coverage_map.py \
    HXB2.bam \
    HXB2.gb \
    HXB2 \
    1 \
    9719 \
    --coverage_scale linear \
    --feature_types CDS,repeat_region,misc_feature,Intron,Exon \
    --feature_colors "CDS=#D3DC7F,repeat_region=green,misc_feature=green,Intron=grey" \
    --output coverage_map.png \