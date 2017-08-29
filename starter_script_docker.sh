pipeline=$1
samples=$2

cmd="luigid --background --pidfile ~/luigi_scheduler_new/luigi_pid --logdir ~/luigi_scheduler_new/ --state-path ~/luigi_scheduler_new/luigi_state &"
echo $cmd
eval $cmd
sleep 4
cd  /srv/qgen/code/qiaseq-singlecell-rna
export LUIGI_CONFIG_PATH=${pipeline}
PYTHONPATH=$PYTHONPATH:"" luigi --module single_cell_rnaseq CombineSamples --samples-cfg ${samples} --worker-wait-interval 20 --workers 16
