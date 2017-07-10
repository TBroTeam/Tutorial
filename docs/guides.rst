Guides
======

Backup data
-----------

Backup your database with these commands::

    docker exec Chado_DB_4_TBro_official pg_dump -U tbro -d chado | xz >tbro_backup_$(date +%F_%T).sql.xz
    docker exec Worker_DB_4_TBro_official pg_dump -U worker -d worker | xz >worker_backup_$(date +%F_%T).sql.xz

You should keep a copy of you blast `.zip` files as well.

Upgrade TBro
------------

If you want to upgrade TBro to a new version (running docker) you have two options.

1. Replace the ``tbro_apache`` container
2. Upgrade in the existing container

.. ATTENTION::
    Those methods will allow you to keep your databases and all of its content.
    However, before you upgrade check for breaking changes in the Changes section of the README.
    And always backup your data. You should do this regularly anyway but especially before performing an upgrade.

Choose option 1 if you did not modify anything inside the ``tbro_apache`` container.
Just execute the following commands::

    docker stop TBro_official
    docker rm TBro_official
    docker pull tbroteam/tbro_apache
    docker run -d --link Chado_DB_4_TBro_official:CHADO --link Worker_FTP_4_TBro_official:WORKERFTP --link Worker_DB_4_TBro_official:WORKER --name "TBro_official" -p 80:80 tbroteam/tbro_apache
    docker exec -i -t TBro_official /home/tbro/build_installation.sh

For option 2 follow these steps::

    docker exec -it TBro_official /bin/bash
    # inside the container
    source ~/.bash_profile
    cd /home
    git clone https://github.com/TBroTeam/TBro.git
    # you now have tbro (old) and TBro (new) in /home
    cd TBro
    # git checkout <branch> # if you want a specific branch instead of master
    cp ../tbro/build.properties .
    phing cli-install
    phing web-install
    phing database-update-modifications


Password protect TBro
---------------------

See the write up by 000generic at https://github.com/TBroTeam/TBro/issues/48


TBro in Docker in Amazon's AWS Lightsail
----------------------------------------

Follow this protocol: https://benchling.com/s/prt-DHSo7HeddC5zM0x7x1Y4
