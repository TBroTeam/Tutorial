Guides
======

Backup data
-----------

Backup your database with these commands::

    docker exec Chado_DB_4_TBro_official pg_dump -U tbro -d chado | xz >tbro_backup_$(date +%F_%T).sql.xz
    docker exec Worker_DB_4_TBro_official pg_dump -U worker -d worker | xz >worker_backup_$(date +%F_%T).sql.xz

You should keep a copy of you blast `.zip` files as well.

Password protect TBro
---------------------

See the write up by 000generic at https://github.com/TBroTeam/TBro/issues/48

