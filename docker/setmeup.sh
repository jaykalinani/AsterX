cd
grep 'module load' .bashrc > /dev/null
if [ $? = 1 ]
then
  echo module load mpi >> .bashrc
fi
source .bashrc
if [ ! -d cactusamrex ]
then
  git clone https://stevenrbrandt@bitbucket.org/eschnett/cactusamrex.git
fi
cd /home/jovyan/cactusamrex
git pull
if [ ! -x GetComponents ]
then
  curl -kLO https://raw.githubusercontent.com/gridaphobe/CRL/ET_2019_03/GetComponents
  chmod a+x GetComponents
  ./GetComponents --parallel azure-pipelines/carpetx.th
fi
mkdir -p /home/jovyan/cactusamrex/Cactus
cd /home/jovyan/cactusamrex/Cactus
echo docker > /home/jovyan/.hostname
if [ ! -r repos/simfactory2/mdb/machines/docker.ini ]
then
  ./simfactory/bin/sim setup-silent
fi
#COPY amrex.cfg ./
./simfactory/bin/sim build -j2 --thornlist ../azure-pipelines/carpetx.th --optionlist /amrex.cfg
#./exe/cactus_sim ./repos/cactusamrex/AMReX/par/amrex.par
