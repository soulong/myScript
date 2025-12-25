################## common terminal ###############
# convert \r problem in Win system
sudo apt install dos2unix
dos2unix xxxx.sh

# release linux cache memory
sudo sh -c 'echo 1 >  /proc/sys/vm/drop_caches'

# check if GPU is okay
python -c "import torch; print(torch.cuda.is_available())"

# count how much files are under current directory
find . -maxdepth 1 -type d -print0 | while IFS= read -r -d '' i ; do
    echo -n $i": " ; 
    (find "$i" -type f | wc -l) ; 
done

# replace all dot with dash of filenames under a directory
for f in *.*; do pre="${f%.*}"; suf="${f##*.}"; mv -i -- "$f" "${pre//./_}.${suf}"; done

# connect to local colab
sudo docker run --rm --runtime=nvidia --gpus all ubuntu nvidia-smi

# rename all file names under a directory
# check first, then remove echo for real running
for i in *target.sorted.bam$*;
do
   echo "mv $i ${i%target.sorted.bam}sorted.bam"
   mv $i ${i%target.sorted.bam}sorted.bam
done

################## WSL ###############
# mount a windows drive to wsl
sudo mkdir /mnt/e
sudo mount -t drvfs E: /mnt/e

# mount a windows share to wsl if already have access to share
# use single quotes otherwise you will need to escape the backslashes
sudo mkdir /mnt/mint
sudo mount -t drvfs '\\10.36.172.157\mint' /mnt/mint

# remove a mounted driver
sudo rm /mnt/driveX

# export/import wsl image
image_name=Ubuntu-24.04
wsl --shutdown && wsl -l -v 
wsl --export $image_name F:/${image_name}.tar
wsl --import $image_name C:/wsl/ D:/${image_name}.tar --version 2
wsl --unregister $image_name # release old image space

# compact wsl image, first shtdown wsl2 by: wsl --shutdown
# in powershell
diskpart
select vdisk file="C:\Users\haohe\AppData\Local\wsl\{4b2b6a94-a636-4186-9794-fb87d75f0b0c}\ext4.vhdx"
attach vdisk readonly
compact vdisk
detach vdisk
exit

/mnt/e/04_LGLab/EED/CutTag_IMR_MEF/MEF

################## NGS ##################
# make sure no backspace in all filepaths
# list files (absolute path)
parameter=$(find . -type f -name "*.filtered.bam" -print0 | xargs -0 realpath | sort | tr '\n' ' ' | sed 's/ $//')
# list files (relative path)
parameter=$(find . -type f -name "*.filtered.bam" -print0 | xargs -0 realpath --relative-to=. | sort | tr '\n' ' ' | sed 's/ $//')
# Convert space-separated string → array using IFS, make sure backspace in file path still works
IFS=' ' read -ra parameterArray <<< "$parameter"
echo "Found ${parameterArray[@]}"

# deeptools heatmap
computeMatrix reference-point -p 8 -R '/mnt/d/Index/hs/v49/genes_protein_coding.bed' \
   -S ${parameterArray[@]} -o genes_protein_coding.gz \
   -a 3000 -b 3000 --skipZeros --smartLabels --quiet
plotHeatmap -m genes_protein_coding.gz -o genes_protein_coding.pdf --colorMap RdBu

# markdup
picard -Xmx48G MarkDuplicates --INPUT $bam --OUTPUT ${prefix}.markdup.bam --REFERENCE_SEQUENCE $fasta --METRICS_FILE ${prefix}_markdup_metric.txt





################## miniforge ##################
# install mamba

# export and import env
mamba env export --no-builds | grep -v "^prefix: " > ngs.yml
mamba env create -f ngs.yml -p ~/miniforge3/envs/ngs

# mamba export all envs
for env in $(mamba env list | cut -d" " -f1); do 
   if [[ ${env:0:1} == "#" ]] ; then continue; fi;
   mamba env export -n $env > ./${env}.yml
done

# mamba remove all installed package in base env
mamba remove `conda list|awk 'NR>3 {print $1}'|tr '\n' ' '`
mamba remove `conda list|awk 'NR>3 {print $1}'|tr '\n' ' '`

    
################## run shinny APP ##################
# run shinyAPP from terminal
cd ~/Documents/GitHub/myScript/R_scripts/shiny_generate_well_positions
nohup R -e "shiny::runApp('shinyApp.R')" &


################## local PC ##################
Anydesk Shulab_Desktop: 1936728718 123456abc-


################## samba share ##################
# start samba server
sudo service smbd restart
sudo ufw allow samba
# list users
sudo pdbedit -L -v 
# add samba user
sudo smbpasswd -a [username]
# change pwd
sudo smbpasswd [username]
# remove samba user
sudo pdbedit -x -u [username]


################## cellprofiler-analyst ##################
# run cellprofiler-analyst
conda activate cellprofiler-analyst
python3 -m CellProfiler-Analyst

# modify sqlite path manually for different OS usage, run below code in SQL supportive software, like sqlite db browser
UPDATE Per_Image SET Image_PathName_ch1="Z:/Data/Project_condensate_reporter/2025-03-30_JY/Images";
UPDATE Per_Image SET Image_PathName_ch2="Z:/Data/Project_condensate_reporter/2025-03-30_JY/Images";
UPDATE Per_Image SET Image_PathName_ch3="Z:/Data/Project_condensate_reporter/2025-03-30_JY/Images";
UPDATE Per_Image SET Image_PathName_ch4="Z:/Data/Project_condensate_reporter/2025-03-30_JY/Images";
UPDATE Per_Image SET Image_ObjectsPathName_mask_cp_masks_cell="Z:/Data/Project_condensate_reporter/2025-03-30_JY/Images";



################## pymol ###############
# pymol reindex
alter (fold_2025_05_19_22_11_atrkd_chk1_sq2_atp_mg_model_0 and chain A),resi=str(int(resi)+2295)
# pymol change chain id
alter (chain B), chain='A'
# pymol align
align mobile, reference
align mobile and chain A, reference and chain B
# remove res
remove sol
remove resn NA,CL
hide hrydrogen
H_add



################## 远程vscode 转发 windows 到 WSL ###############

# 以管理员身份运行
$wslIp = (wsl hostname -I).Trim()
if ($wslIp -eq "") {
    Write-Host "WSL 未运行或无法获取 IP"
    exit 1
}

# 删除旧规则
netsh interface portproxy delete v4tov4 listenport=2222 listenaddress=0.0.0.0

# 添加新规则
netsh interface portproxy add v4tov4 listenport=2222 listenaddress=0.0.0.0 connectport=2222 connectaddress=$wslIp

# 确保防火墙放行（可加判断）
if (-not (Get-NetFirewallRule -DisplayName "WSL2 SSH" -ErrorAction SilentlyContinue)) {
    New-NetFirewallRule -DisplayName "WSL2 SSH" -Direction Inbound -Protocol TCP -LocalPort 2222 -Action Allow
}

Write-Host "✅ 已设置端口转发: Windows:2222 → WSL($wslIp):2222"
