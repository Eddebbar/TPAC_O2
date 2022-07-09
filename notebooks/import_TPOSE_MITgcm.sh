#!/bin/bash

set -e

local_dir=/glade/scratch/yeddebba/TPOSE/

if ! [ -e ${local_dir} ]
then
    mkdir ${local_dir}
else
    echo ${local_dir} "exists!"     
fi

remote_dir=http://www.ecco.ucsd.edu/DATA/TROPAC/tpose6
wget -r -l1 --no-parent -A ".nc" ${remote_dir} -P ${local_dir}

remote_dir=http://www.ecco.ucsd.edu/DATA/TROPAC/bgc/monthly
wget -r -l1 --no-parent -A ".nc" ${remote_dir} -P ${local_dir}

remote_dir=http://www.ecco.ucsd.edu/DATA/TROPAC/bgc
wget -r -l1 --no-parent -A ".nc" ${remote_dir} -P ${local_dir}
