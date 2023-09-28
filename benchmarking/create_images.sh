#!/bin/bash
#SBATCH --cpus-per-task=18
#SBATCH --job-name=create_images
#SBATCH --mem=250G
#SBATCH --partition=slim16
#SBATCH --time=12-00:00:00
#SBATCH -o slurm-out-create_images-%j.out

RADICL_INDEX=radicl-mouse/index.50000.smoother_index
HIC_INDEX=two_replicate_index.smoother_index

EXPORT_PARAMS="-f png -s 2500"

create_index() {
    biosmoother init ${HIC_INDEX} organisms/t_brucei/genome/HGAP3_Tb427v10.sizes organisms/t_brucei/genome/HGAP3_Tb427v10.gff -d 10000

    zcat organisms/t_brucei/data/merged.1.pairs.gz | biosmoother repl ${HIC_INDEX} - WT -g a
    zcat organisms/t_brucei/data/SRR7721307.1.pairs.gz | biosmoother repl ${HIC_INDEX} - H3.V-/- -g b

    module load ngs/bwa/0.7.16
    bwa mem -t 18 organisms/t_brucei/genome/HGAP3_Tb427v10.fna organisms/t_brucei/fasta/disabled/SRR7988175_R1.fq.gz organisms/t_brucei/fasta/disabled/SRR7988175_R2.fq.gz | awk '!/^#|^@/ {print $1, $3, $4, $2 % 16 == 0 ? "+" :"-", $5, $12}' OFS="\t" > organisms/t_brucei/data/disabled/SRR7988175.tsv
    biosmoother track -g col ${HIC_INDEX} organisms/t_brucei/data/disabled/SRR7988175.tsv chip-seq
}
# create_index

default_hi_c () {
    libbiosmoother reset $HIC_INDEX
    libbiosmoother set $HIC_INDEX settings.filters.symmetry mirror

    libbiosmoother set $HIC_INDEX settings.interface.show_hide.axis false
    libbiosmoother set $HIC_INDEX settings.interface.show_hide.regs false
    libbiosmoother set $HIC_INDEX settings.interface.show_hide.raw false
    libbiosmoother set $HIC_INDEX settings.interface.show_hide.annotation false

    libbiosmoother set $HIC_INDEX settings.filters.show_contig_smaller_than_bin true
    libbiosmoother set $HIC_INDEX settings.filters.cut_off_bin fit_chrom_larger

    libbiosmoother set $HIC_INDEX settings.normalization.normalize_by ice
    libbiosmoother set $HIC_INDEX settings.normalization.ice_sparse_slice_filter.val 0
    libbiosmoother set $HIC_INDEX settings.export.print_region false
}

show_chr_8 () {
    libbiosmoother set $HIC_INDEX area.x_end 2341
    libbiosmoother set $HIC_INDEX area.x_start 2121
    libbiosmoother set $HIC_INDEX area.y_end 2341
    libbiosmoother set $HIC_INDEX area.y_start 2121
}

fig1 () {
    # Fig. 1 B left panel, top: whole dataset
    echo "Fig. 1 B left panel, top: whole dataset"
    default_hi_c
    show_chr_8
    libbiosmoother export $HIC_INDEX -p pics/fig1b_left_top $EXPORT_PARAMS

    # Fig. 1 B left panel, middle: mapping quality
    echo "Fig. 1 B left panel, middle: mapping quality"
    default_hi_c
    show_chr_8
    libbiosmoother set $HIC_INDEX settings.filters.mapping_q.val_min 2
    libbiosmoother export $HIC_INDEX -p pics/fig1b_left_middle $EXPORT_PARAMS

    # Fig. 1 B left panel, bottom: annotation filter
    echo "Fig. 1 B left panel, bottom: annotation filter"
    default_hi_c
    show_chr_8
    libbiosmoother set $HIC_INDEX annotation.filter_absent_x "gene"
    libbiosmoother export $HIC_INDEX -p pics/fig1b_left_bottom $EXPORT_PARAMS


    # Fig. 1 B center panel, top: 10 kbp
    echo "Fig. 1 B center panel, top: 10 kbp"
    default_hi_c
    show_chr_8
    libbiosmoother set $HIC_INDEX settings.interface.fixed_bin_size true
    libbiosmoother set $HIC_INDEX settings.interface.fixed_bin_size_x.val 1
    libbiosmoother set $HIC_INDEX settings.interface.fixed_bin_size_y.val 1
    libbiosmoother export $HIC_INDEX -p pics/fig1b_center_top $EXPORT_PARAMS

    # Fig. 1 B center panel, middle: 50 kbp
    echo "Fig. 1 B center panel, middle: 50 kbp"
    default_hi_c
    show_chr_8
    libbiosmoother set $HIC_INDEX settings.interface.fixed_bin_size true
    libbiosmoother set $HIC_INDEX settings.interface.fixed_bin_size_x.val 5
    libbiosmoother set $HIC_INDEX settings.interface.fixed_bin_size_y.val 5
    libbiosmoother export $HIC_INDEX -p pics/fig1b_center_middle $EXPORT_PARAMS

    # Fig. 1 B center panel, bottom: 200 kbp
    echo "Fig. 1 B center panel, bottom: 200 kbp"
    default_hi_c
    show_chr_8
    libbiosmoother set $HIC_INDEX settings.interface.fixed_bin_size true
    libbiosmoother set $HIC_INDEX settings.interface.fixed_bin_size_x.val 20
    libbiosmoother set $HIC_INDEX settings.interface.fixed_bin_size_y.val 20
    libbiosmoother export $HIC_INDEX -p pics/fig1b_center_bottom $EXPORT_PARAMS


    # Fig. 1 B right panel, top: 10 kbp
    echo "Fig. 1 B right panel, top: 10 kbp"
    default_hi_c
    show_chr_8
    libbiosmoother set $HIC_INDEX settings.normalization.normalize_by dont
    libbiosmoother export $HIC_INDEX -p pics/fig1b_right_top $EXPORT_PARAMS

    # Fig. 1 B right panel, middle: 50 kbp
    echo "Fig. 1 B right panel, middle: 50 kbp"
    default_hi_c
    show_chr_8
    libbiosmoother export $HIC_INDEX -p pics/fig1b_right_middle $EXPORT_PARAMS

    # Fig. 1 B right panel, bottom: 200 kbp
    echo "Fig. 1 B right panel, bottom: 200 kbp"
    default_hi_c
    show_chr_8
    libbiosmoother set $HIC_INDEX settings.normalization.ddd true
    libbiosmoother export $HIC_INDEX -p pics/fig1b_right_bottom $EXPORT_PARAMS
}
# fig1




default_radicl () {
    libbiosmoother reset $RADICL_INDEX

    libbiosmoother set $RADICL_INDEX settings.interface.show_hide.axis false
    libbiosmoother set $RADICL_INDEX settings.interface.show_hide.regs false
    libbiosmoother set $RADICL_INDEX settings.interface.show_hide.raw false
    libbiosmoother set $RADICL_INDEX settings.interface.show_hide.annotation false

    libbiosmoother set $RADICL_INDEX settings.filters.show_contig_smaller_than_bin true
    libbiosmoother set $RADICL_INDEX settings.filters.cut_off_bin cover_multiple
    libbiosmoother set $RADICL_INDEX settings.export.print_region false
}

fig2 () {
    # Fig. 2 A top: any mappging quality
    echo "Fig. 2 A top: any mappging quality"
    default_radicl
    libbiosmoother set $RADICL_INDEX area.x_end 3000
    libbiosmoother set $RADICL_INDEX area.x_start 180
    libbiosmoother set $RADICL_INDEX area.y_end 3000
    libbiosmoother set $RADICL_INDEX area.y_start 180
    libbiosmoother export $RADICL_INDEX -p pics/fig2a_top $EXPORT_PARAMS

    # Fig. 2 A bottom: high mappging quality
    echo "Fig. 2 A bottom: high mappging quality"
    default_radicl
    libbiosmoother set $RADICL_INDEX area.x_end 3000
    libbiosmoother set $RADICL_INDEX area.x_start 180
    libbiosmoother set $RADICL_INDEX area.y_end 3000
    libbiosmoother set $RADICL_INDEX area.y_start 180
    libbiosmoother set $RADICL_INDEX settings.filters.mapping_q.val_min 2
    libbiosmoother export $RADICL_INDEX -p pics/fig2a_bottom $EXPORT_PARAMS

    # Fig. 2 B top: all multimappers
    echo "Fig. 2 B top: all multimappers"
    default_hi_c
    libbiosmoother set $HIC_INDEX area.x_end 722
    libbiosmoother set $HIC_INDEX area.x_start 577
    libbiosmoother set $HIC_INDEX area.y_end 722
    libbiosmoother set $HIC_INDEX area.y_start 577
    libbiosmoother set $HIC_INDEX settings.filters.ambiguous_mapping overlaps
    libbiosmoother export $HIC_INDEX -p pics/fig2b_top $EXPORT_PARAMS

    echo "Fig. 2 B bottom: enclosed multimappers"
    default_hi_c
    libbiosmoother set $HIC_INDEX area.x_end 722
    libbiosmoother set $HIC_INDEX area.x_start 577
    libbiosmoother set $HIC_INDEX area.y_end 722
    libbiosmoother set $HIC_INDEX area.y_start 577
    libbiosmoother export $HIC_INDEX -p pics/fig2b_bottom $EXPORT_PARAMS

    # Fig. 2 C top: any position
    echo "Fig. 2 C top: any position"
    default_radicl
    libbiosmoother set $RADICL_INDEX area.x_end 2540
    libbiosmoother set $RADICL_INDEX area.x_start 1364
    libbiosmoother set $RADICL_INDEX area.y_end 2540
    libbiosmoother set $RADICL_INDEX area.y_start 1364
    libbiosmoother export $RADICL_INDEX -p pics/fig2c_top $EXPORT_PARAMS

    # Fig. 2 C middle: overlapping gene
    echo "Fig. 2 C top: overlapping gene"
    default_radicl
    libbiosmoother set $RADICL_INDEX area.x_end 2540
    libbiosmoother set $RADICL_INDEX area.x_start 1364
    libbiosmoother set $RADICL_INDEX area.y_end 2540
    libbiosmoother set $RADICL_INDEX area.y_start 1364
    libbiosmoother set $RADICL_INDEX annotation.filter_absent_x gene
    libbiosmoother export $RADICL_INDEX -p pics/fig2c_middle $EXPORT_PARAMS

    # Fig. 2 C bottom: gene coordinates
    echo "Fig. 2 C bottom: gene coordinates"
    default_radicl
    libbiosmoother set $RADICL_INDEX area.x_end 1097
    libbiosmoother set $RADICL_INDEX area.x_start 642
    libbiosmoother set $RADICL_INDEX area.y_end 2540
    libbiosmoother set $RADICL_INDEX area.y_start 1364
    libbiosmoother set $RADICL_INDEX settings.filters.anno_coords_col true
    libbiosmoother export $RADICL_INDEX -p pics/fig2c_bottom $EXPORT_PARAMS

    # ploidy coords
    echo "Fig. 2 D top: uncorrected"
    default_hi_c
    biosmoother ploidy ${HIC_INDEX} organisms/t_brucei/genome/HGAP3_Tb427v10.ploidy
    libbiosmoother set $HIC_INDEX area.x_end 1206
    libbiosmoother set $HIC_INDEX area.x_start 756
    libbiosmoother set $HIC_INDEX area.y_end 1206
    libbiosmoother set $HIC_INDEX area.y_start 756
    libbiosmoother set $HIC_INDEX settings.normalization.log_base.val 10
    libbiosmoother set $HIC_INDEX settings.normalization.ploidy_correct false
    libbiosmoother export $HIC_INDEX -p pics/fig2d_top $EXPORT_PARAMS

    echo "Fig. 2 D bottom: ploidy corrected"
    default_hi_c
    biosmoother ploidy ${HIC_INDEX} organisms/t_brucei/genome/HGAP3_Tb427v10.ploidy
    libbiosmoother set $HIC_INDEX area.x_end 1206
    libbiosmoother set $HIC_INDEX area.x_start 756
    libbiosmoother set $HIC_INDEX area.y_end 1206
    libbiosmoother set $HIC_INDEX area.y_start 756
    libbiosmoother set $HIC_INDEX settings.normalization.log_base.val 10
    libbiosmoother export $HIC_INDEX -p pics/fig2d_bottom $EXPORT_PARAMS

    # comparison
    echo "Fig. 2 E: comparison"
    default_hi_c
    libbiosmoother set $HIC_INDEX area.x_end 1233
    libbiosmoother set $HIC_INDEX area.x_start 859
    libbiosmoother set $HIC_INDEX area.y_end 1233
    libbiosmoother set $HIC_INDEX area.y_start 859
    libbiosmoother set $HIC_INDEX settings.replicates.between_group sub
    libbiosmoother set $HIC_INDEX settings.normalization.color_range.val_min -1
    libbiosmoother set $HIC_INDEX settings.normalization.scale abs
    libbiosmoother set $HIC_INDEX settings.normalization.log_base.val 15
    libbiosmoother export $HIC_INDEX -p pics/fig2e $EXPORT_PARAMS

    # v4c 1
    echo "Fig. 2 F: V4C"
    default_hi_c
    libbiosmoother set $HIC_INDEX area.x_end 1960
    libbiosmoother set $HIC_INDEX area.x_start 1732
    libbiosmoother set $HIC_INDEX area.y_end 4065
    libbiosmoother set $HIC_INDEX area.y_start 3837
    libbiosmoother set $HIC_INDEX settings.interface.show_hide.raw true
    libbiosmoother set $HIC_INDEX settings.interface.show_hide.annotation true
    libbiosmoother set $HIC_INDEX annotation.visible_x ""
    libbiosmoother set $HIC_INDEX annotation.visible_y ""
    libbiosmoother set $HIC_INDEX settings.interface.max_num_bins.val 10
    libbiosmoother set $HIC_INDEX settings.interface.v4c.do_col true
    libbiosmoother set $HIC_INDEX settings.interface.v4c.col_from 3941
    libbiosmoother set $HIC_INDEX settings.interface.v4c.col_to 3956
    libbiosmoother set $HIC_INDEX settings.interface.raw_size.val 2500
    libbiosmoother set $HIC_INDEX settings.interface.anno_size.val 100
    libbiosmoother set $HIC_INDEX settings.interface.center_tracks_on_bins true
    libbiosmoother set $HIC_INDEX settings.interface.connect_tracks_over_contig_borders true
    libbiosmoother set $HIC_INDEX settings.normalization.log_base.val 3
    libbiosmoother set $HIC_INDEX settings.export.secondary_stroke_width.val 50
    libbiosmoother set $HIC_INDEX settings.export.white_background false
    libbiosmoother set $HIC_INDEX settings.export.secondary_x_range.val_min 0
    libbiosmoother set $HIC_INDEX settings.export.secondary_x_range.val_max 0.007
    libbiosmoother set $HIC_INDEX coverage.in_column ""
    biosmoother export $HIC_INDEX -p pics/fig2f $EXPORT_PARAMS

    # v4c 2
    echo "Fig. 2 F: Extra V4C"
    default_hi_c
    libbiosmoother set $HIC_INDEX area.x_end 1960
    libbiosmoother set $HIC_INDEX area.x_start 1732
    libbiosmoother set $HIC_INDEX area.y_end 4065
    libbiosmoother set $HIC_INDEX area.y_start 3837
    libbiosmoother set $HIC_INDEX settings.interface.show_hide.raw true
    libbiosmoother set $HIC_INDEX settings.interface.show_hide.annotation true
    libbiosmoother set $HIC_INDEX annotation.visible_x ""
    libbiosmoother set $HIC_INDEX annotation.visible_y ""
    libbiosmoother set $HIC_INDEX settings.interface.max_num_bins.val 10
    libbiosmoother set $HIC_INDEX settings.interface.v4c.do_col true
    libbiosmoother set $HIC_INDEX settings.interface.v4c.col_from 3881
    libbiosmoother set $HIC_INDEX settings.interface.v4c.col_to 3896
    libbiosmoother set $HIC_INDEX settings.interface.raw_size.val 2500
    libbiosmoother set $HIC_INDEX settings.interface.anno_size.val 100
    libbiosmoother set $HIC_INDEX settings.interface.center_tracks_on_bins true
    libbiosmoother set $HIC_INDEX settings.interface.connect_tracks_over_contig_borders true
    libbiosmoother set $HIC_INDEX settings.normalization.log_base.val 3
    libbiosmoother set $HIC_INDEX settings.export.secondary_stroke_width.val 50
    libbiosmoother set $HIC_INDEX settings.export.white_background false
    libbiosmoother set $HIC_INDEX settings.export.secondary_x_range.val_min 0
    libbiosmoother set $HIC_INDEX settings.export.secondary_x_range.val_max 0.007
    libbiosmoother set $HIC_INDEX coverage.in_column ""
    biosmoother export $HIC_INDEX -p pics/fig2f_extra $EXPORT_PARAMS

    echo "Fig. 2 G: overlays"
    default_hi_c
    libbiosmoother set $HIC_INDEX settings.interface.show_hide.raw true
    libbiosmoother set $HIC_INDEX settings.interface.show_hide.annotation true
    libbiosmoother set $HIC_INDEX area.x_end 1911 # chr7_core 2,270 kbp # chr7_5A is 860kbp long
    libbiosmoother set $HIC_INDEX area.x_start 1842 # chr7_5A 1,580 kbp
    libbiosmoother set $HIC_INDEX area.y_end 1911
    libbiosmoother set $HIC_INDEX area.y_start 1842
    libbiosmoother set $HIC_INDEX settings.interface.raw_size.val 2000
    libbiosmoother set $HIC_INDEX settings.interface.anno_size.val 500
    libbiosmoother set $HIC_INDEX annotation.visible_x "gene"
    libbiosmoother set $HIC_INDEX annotation.visible_y "gene"
    libbiosmoother set $HIC_INDEX settings.export.secondary_stroke_width.val 50
    libbiosmoother set $HIC_INDEX settings.export.white_background false
    libbiosmoother set $HIC_INDEX settings.export.secondary_x_range.val_min 0
    libbiosmoother set $HIC_INDEX settings.export.secondary_x_range.val_max 15500
    biosmoother export $HIC_INDEX -p pics/fig2g $EXPORT_PARAMS

    # Extra pictures: assoc. slices norm
    echo "Extra: assoc. slices"
    default_radicl
    libbiosmoother set $RADICL_INDEX area.x_end 3000
    libbiosmoother set $RADICL_INDEX area.x_start 180
    libbiosmoother set $RADICL_INDEX area.y_end 3000
    libbiosmoother set $RADICL_INDEX area.y_start 180
    libbiosmoother set $RADICL_INDEX settings.normalization.normalize_by grid-seq
    libbiosmoother export $RADICL_INDEX -p pics/extra_1_assoc_slices $EXPORT_PARAMS

    # Extra pictures: binom. test
    echo "Extra: binom. test"
    default_radicl
    libbiosmoother set $RADICL_INDEX area.x_end 3000
    libbiosmoother set $RADICL_INDEX area.x_start 180
    libbiosmoother set $RADICL_INDEX area.y_end 3000
    libbiosmoother set $RADICL_INDEX area.y_start 180
    libbiosmoother set $RADICL_INDEX settings.normalization.normalize_by radicl-seq
    libbiosmoother export $RADICL_INDEX -p pics/extra_1_binom_test $EXPORT_PARAMS
}
fig2


echo "done"
