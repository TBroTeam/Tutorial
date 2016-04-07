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


Preload demo data
-----------------

