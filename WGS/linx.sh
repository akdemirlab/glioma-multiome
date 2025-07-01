#BSUB -W 240:00
#BSUB -q e40long
#BSUB -n 12
#BSUB -M 475
#BSUB -R rusage[mem=475]
#BSUB -J $<tumor>
#BSUB -cwd $<MSUB_OUTPUT_DIR>
#BSUB -o $<MSUB_OUTPUT_DIR>/$<tumor>.out

. ~/.bashrc
module load openjdk/11.0.5-10
module load R/4.2.1
module load bwa/0.7.17
module load samtools/1.15
module load perl/5.28.1
export PATH=/rsrch4/home/neurosurgery/akdemirlab/pipelines/circos-0.69-9/bin:$PATH

script=/rsrch4/scratch/neurosurgery/akdemirlab/catalyst_wgs/scripts/sv/linx_run_pipeline.sh
pipeline=/rsrch4/home/neurosurgery/akdemirlab/pipelines/linx/pipeline
resources=/rsrch4/home/neurosurgery/akdemirlab/pipelines/linx/resources
tools=/rsrch4/home/neurosurgery/akdemirlab/pipelines/linx/tools

data_dir=/rsrch4/scratch/neurosurgery/akdemirlab/catalyst_wgs/data/bwa_bam_new01
result_dir=/rsrch4/scratch/neurosurgery/akdemirlab/catalyst_wgs/results_new01/linx

cd /rsrch4/scratch/neurosurgery/akdemirlab/catalyst_wgs/results_new01/
$script ${pipeline}/scripts ${data_dir} ${resources} ${tools} '$<tumor>,$<normal>' V38 WGS 12 256 $result_dir
