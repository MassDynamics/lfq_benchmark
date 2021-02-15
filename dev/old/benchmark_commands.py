# iPRG2015
from benchmarkers import IPRG2015Benchmarker
from confusion_matrix_calculator import ConfusionMatrixCalculator   
experiment_home = "../../../experiments/IPRG2015_NewWorkflow" #before tfd-ml
bm = IPRG2015Benchmarker(experiment_home)
bm.run()
print(ConfusionMatrixCalculator().run(bm))
experiment_home = "/mnt/d/Dropbox/MassDynamics_local/experiments/IPRG2015_REPRISAL_2" #not using razor peptides for quant. 
bm = IPRG2015Benchmarker(experiment_home)
bm.run()
print(ConfusionMatrixCalculator().run(bm))
experiment_home = "/mnt/d/Dropbox/MassDynamics_local/experiments/IPRG2015_REPRISAL" #razor peptides for quant
bm = IPRG2015Benchmarker(experiment_home)
bm.run()
print(ConfusionMatrixCalculator().run(bm))
experiment_home = "/mnt/d/Dropbox/MassDynamics_local/experiments/IPRG2015_Reprisal_old_tfd_ml" #new PI, oldtfdML
bm = IPRG2015Benchmarker(experiment_home)
bm.run()
print(ConfusionMatrixCalculator().run(bm))
experiment_home = "/mnt/d/Dropbox/MassDynamics_local/experiments/IPRG2015_REPRISAL" #razor peptides for quant
bm = IPRG2015Benchmarker(experiment_home)
bm.run()
print(ConfusionMatrixCalculator().run(bm))

# Proteomics Human vs Yeast 

from benchmarkers import ProteomicMixtureBenchmark

experiment_home = "/mnt/d/Dropbox/MassDynamics_local/experiments/PC3_Yeast_Original"
ProteomicMixtureBenchmark(experiment_home).run()

experiment_home = "/mnt/d/Dropbox/MassDynamics_local/experiments/PC3_Yeast_REPRISAL"
ProteomicMixtureBenchmark(experiment_home).run()


# MQ, BYO, PD

from benchmarkers import IPRG2015Benchmarker
from confusion_matrix_calculator import ConfusionMatrixCalculator

experiment_home = "/mnt/d/Dropbox/MassDynamics_local/experiments/IPRG2015_REPRISAL" #razor peptides for quant
bm = IPRG2015Benchmarker(experiment_home)
bm.run()
print(ConfusionMatrixCalculator().run(bm))

mq_home = "/mnt/d/Dropbox/MassDynamics_local/experiments/IPRG2015_MQ_BYO/"
bm = IPRG2015Benchmarker(mq_home, "BYO")
bm.run()

mq_home = "/mnt/d/Dropbox/MassDynamics_local/experiments/IPRG2015_MQ_Perseus/"
bm = IPRG2015Benchmarker(mq_home, "Perseus")
bm.run()
print(ConfusionMatrixCalculator().run(bm))

experiment_home = "/mnt/d/Dropbox/MassDynamics_local/experiments/IPRG2015_PD"
bm = IPRG2015Benchmarker_PD(experiment_home)
bm.run()
print(ConfusionMatrixCalculator().run(bm))

# new tfd thresholds

from benchmarkers import IPRG2015Benchmarker     
from confusion_matrix_calculator import ConfusionMatrixCalculator

experiment_home = "/mnt/d/Dropbox/MassDynamics_local/experiments/tfd_tests/prob_05" 
bm = IPRG2015Benchmarker(experiment_home)
bm.run()
print(ConfusionMatrixCalculator().run(bm))

experiment_home = "/mnt/d/Dropbox/MassDynamics_local/experiments/tfd_tests/prob_07" 
bm = IPRG2015Benchmarker(experiment_home)
bm.run()

experiment_home = "/mnt/d/Dropbox/MassDynamics_local/experiments/tfd_tests/prob_062" 
bm = IPRG2015Benchmarker(experiment_home)
bm.run()
print(ConfusionMatrixCalculator().run(bm))

# UPS Benchmarking 

from benchmarkers import IPRG2015Benchmarker , UPSBenchmarker  
from confusion_matrix_calculator import ConfusionMatrixCalculator

experiment_home = "/mnt/d/Dropbox/MassDynamics_local/experiments/MaxLFQ_USP_current_workflow"
bm = UPSBenchmarker(experiment_home)
bm.run()
print(ConfusionMatrixCalculator().run(bm))

experiment_home = "/mnt/d/Dropbox/MassDynamics_local/experiments/ups/clf_prediction"
bm = UPSBenchmarker(experiment_home)
bm.run()
print(ConfusionMatrixCalculator().run(bm))

experiment_home = "/mnt/d/Dropbox/MassDynamics_local/experiments/ups/clf_prob_05"
bm = UPSBenchmarker(experiment_home)
bm.run()
print(ConfusionMatrixCalculator().run(bm))

experiment_home = "/mnt/d/Dropbox/MassDynamics_local/experiments/ups/clf_prob_07"
bm = UPSBenchmarker(experiment_home)
bm.run()
print(ConfusionMatrixCalculator().run(bm))

mq_home = "/mnt/d/Dropbox/MassDynamics_local/experiments/ups_byo/"
bm = UPSBenchmarker(mq_home, mode = "BYO")
bm.run()
print(ConfusionMatrixCalculator().run(bm))

mq_home = "/mnt/d/Dropbox/MassDynamics_local/experiments/ups_perseus/"
bm = UPSBenchmarker(mq_home, mode = "Perseus")
bm.run()
print(ConfusionMatrixCalculator().run(bm))