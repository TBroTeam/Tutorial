Installation
============

Installation of TBro is straightforward and easy using preconfigured `docker <https://www.docker.com/>`_ containers.
See `docker documentation <https://docs.docker.com/engine/installation/>`_ on how to install docker on your machine.

After installation of docker execute the following commands to pull the docker images:
 
::                                                                                                                                          
                                                                                                                                            
        docker pull greatfireball/generic_postgresql_db
        docker pull tbroteam/generic_chado_db_reload
        docker pull tbroteam/tbro_worker_ftp
        docker pull tbroteam/tbro_worker
        docker pull tbroteam/tbro_apache
        docker pull tbroteam/tbro_demo

Now start the Chado database container and install the schema:

::

    docker run -d -e DB_NAME=chado -e DB_USER=tbro -e DB_PW=tbro --name "Chado_DB_4_TBro_official" greatfireball/generic_postgresql_db
    sleep 60
    docker run --rm -i -t --link Chado_DB_4_TBro_official:CHADO --name "Chado_DB_4_TBro_load_official" tbroteam/generic_chado_db_reload

Now start the database container for the BLAST worker:

::

    docker run -d -e DB_NAME=worker -e DB_USER=worker -e DB_PW=worker --name "Worker_DB_4_TBro_official" greatfireball/generic_postgresql_db

Start an ftp server to host the BLAST databases:

::

    docker run -d --name "Worker_FTP_4_TBro_official" -e FTP_USER="tbro" -e FTP_PW="ftp" tbroteam/tbro_worker_ftp

Start a worker to execute the BLAST jobs:

::

    docker run -d --link Worker_DB_4_TBro_official:WORKER --link Worker_FTP_4_TBro_official:WORKERFTP --name "TBro_Worker_official" tbroteam/tbro_worker
    docker exec -i -t TBro_Worker_official /home/tbro/worker_build_installation.sh 2> run_worker_build_installation.err > run_worker_build_installation.log

Finally start and install the main TBro container:

::

    docker run -d --link Chado_DB_4_TBro_official:CHADO --link Worker_FTP_4_TBro_official:WORKERFTP --link Worker_DB_4_TBro_official:WORKER --name "TBro_official" -p 80:80 tbroteam/tbro_apache
    docker exec -i -t TBro_official /home/tbro/build_installation.sh 2> tbro_build_installation.err > tbro_build_installation.log

You can now access the TBro web interface by pointing your browser to http://localhost 
However there is no data loaded, yet.
To load data you can either perform the automatic demo installation (see next section), follow the step-by-step tutorial (next chapter) or load your own data.

Preload demo data
-----------------

Run the following command to fill your TBro instance with demo data from *Cannabis sativa*.

::

    docker run --rm -i -t --link Worker_DB_4_TBro_official:WORKER --link Worker_FTP_4_TBro_official:WORKERFTP --link Chado_DB_4_TBro_official:CHADO --name "TBro_Demo_official" tbroteam/tbro_demo

Congratulations, you have a full-featured TBro instance up and running.

