# CLIPA
Patient-derived Cancer Cell Line Interactive Pharmacogenomic Analyses

Welcome to view source code for CLIPA!

!!! ALL necessary files and scripts is open for every users. However, due to file size limitation, the dataset( mainly *.RData ) are not here. One get the urls of all RDatas in urls.txt.

CLIPA dependents on php and Rscript. 

It could running on any web server, like apache and IIS.

Steps to get CLIPA running on Debian 10 for example:

1. Install apache2:

   apt install apache2
   
2. Install dependencys:

   apt -y install php-cli php-fpm php-json php-pdo php-mysql php-zip php-gd php-mbstring php-curl php-xml php-pear php-bcmath libapache2-mod-php libxml2-dev libcurl4-openssl-dev
   
3. Start apache2:

   service apache2 start
   
4. Enable php7.3:
   
   a2enmod php7.3
   
   systemctl restart apache2

5. Install latest R:

   NOTE: 'apt -y install r-base' will not install the latest R. One can follow the steps below to install the latest R:
      
      1. vim /etc/apt/sources.list, and add the following repository:
       
         deb http://cloud.r-project.org/bin/linux/debian buster-cran40/
         
      2. Update your repository list:
      
         apt update
 
         This step will return a warning. Find and copy the key after 'NO_PUBKEY' promote, like FCAE2A0E115C3D8A.
         
      3. Add key:
         
         apt-key adv --keyserver keyserver.ubuntu.com --recv-keys FCAE2A0E115C3D8A
         
      4. Install R
      
         apt install -t buster-cran40 r-base
         
      5. Input 'R' to enter R console and install the following packages:
      
         rjson, GSVA, Dplyr, Impute, Cluster, reshape2, tidyr, stringr, caret, randomForest, e1071, kknn, rpart, gbm, adabag, glmnet, nnet, pROC
         
6. Put source codes into /var/www/html/
         
         
Hope you enjoy it!
         
