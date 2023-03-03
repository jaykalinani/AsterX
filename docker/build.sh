PROCS=$(lscpu | grep CPU.s.: | head -1 | cut -d: -f2)
if [ "$USE_PROCS" = "" ]
then
  USE_PROCS=$(($PROCS/2))
fi
cd
if [ ! -d cactusamrex ]
then
  git clone https://bitbucket.org/eschnett/cactusamrex.git
fi
cd /home/jovyan/cactusamrex
git pull
if [ ! -x GetComponents ]
then
  wget https://raw.githubusercontent.com/gridaphobe/CRL/ET_2019_03/GetComponents
  chmod a+x GetComponents
  ./GetComponents --parallel azure-pipelines/carpetx.th
fi
cd /home/jovyan/cactusamrex/Cactus
echo docker > /home/jovyan/.hostname
if [ ! -r repos/simfactory2/mdb/machines/docker.ini ]
then
  ./simfactory/bin/sim setup-silent
fi
./simfactory/bin/sim build -j ${USE_PROCS} --thornlist ../azure-pipelines/carpetx.th --optionlist ../azure-pipelines/ubuntu.cfg 2>&1 | tee make.out
#./exe/cactus_sim ./repos/cactusamrex/AMReX/par/amrex.par
