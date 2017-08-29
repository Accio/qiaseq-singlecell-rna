pipeline-cfg=$1
samples-cfg=$2

luigid --background --pidfile ~/luigi_scheduler_new/luigi_pid --logdir ~/luigi_scheduler_new/ --state-path ~/luigi_scheduler_new/luigi_state &
cd  /srv/qgen/code/
git clone https://containerUser:okazaki1@github.com/qiaseq/qiaseq-singlecell-rna.git
cd qiaseq-singlecell-rna
export LUIGI_CONFIG_PATH=${pipeline-cfg}
PYTHONPATH=$PYTHONPATH:"" luigi --module CombineSamples --samples-cfg ${samples-cfg} --worker-wait-interval 20 --workers 16
