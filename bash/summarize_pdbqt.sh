#! /bin/bash
# env: unidock
# need ringtail installed


working_dir="/home/hao/docking_tools/VS_HPK1"

docking_pdbqt_output_dir="vina_result_mcule_ue1"
recepor_file="7M0M.pdbqt"
output_db="vina_result_mcule_ue1.db"


cd $working_dir

rt_process_vs write \
    --file_path $docking_pdbqt_output_dir \
    --output_db $output_db \
    --docking_mode vina \
    --receptor $recepor_file \
    --save_receptor \
    --max_poses 1 \
    --print_summary
    # --input_db xxx --append_results
    # --overwrite

# rt_process_vs read --input_db $output_db --score_percentile 1 --bookmark_name score_p1

# rt_process_vs read --input_db $output_db --bookmark_name score_p5 --export_sdf_path score_p5 --pymol