#!/usr/bin/env sh

################################################################################
# File: checkout_data_files.sh
#
# The purpose of this script is to obtain lookup tables used by the WRF physics
#   packages.  At present, the only method for acquiring these tables is through
#   the iAMAS-Dev GitHub repository using either git, svn, or curl.
#
# If none of the methods used in this script are successful in acquiring the 
#   tables, please attempt to manually download the files from the iAMAS-Data 
#   repository at https://github.com/iAMAS-Dev/iAMAS-Data/.  All *.TBL and *DATA* 
#   files, as well as the COMPATIBILITY file, should be copied into 
#   a subdirectory named src/core_atmosphere/physics/physics_wrf/files before 
#   continuing the build process.  In general, one should obtain the lookup
#   tables from a tag in the iAMAS-Dev repository whose name matches the version
#   of the iAMAS-Atmosphere code; e.g., for iAMAS-Atmosphere v7.0, one should
#   use the tables from the v7.0 tag in the iAMAS-Data repository.
#
# If all else fails, please contact the iAMAS-Atmosphere developers through
#   the iAMAS-Atmosphere support forum at http://forum.mmm.ucar.edu/.
#
################################################################################


iamas_vers="7.0"

github_org="iAMAS-Dev"   # GitHub organization where the iAMAS-Data repository is found.
                        # For physics development, it can be helpful for a developer
                        # to obtain tables from their own fork of the iAMAS-Data repository.

#
# Return 1 if the "iamas_vers" string is found in the physics table COMPATIBILITY
# file, and 0 otherwise
#
check_compatibility() {
   for ver in `cat physics_wrf/files/COMPATIBILITY | grep -v "#"`; do
      if [ "${ver}" = "${iamas_vers}" ]; then
         return 1
      fi
   done
   return 0
}


#
# See whether we already have compatible physics tables
#
if [ -s physics_wrf/files/COMPATIBILITY ]; then

   check_compatibility
   if [ $? -eq 1 ]; then
      echo "*** Compatible versions of WRF physics tables appear to already exist; no need to obtain them again ***"
      exit 0
   else
      echo "*** Existing WRF physics tables appear to be incompatible with iAMAS v${iamas_vers}; attempting to download compatible tables ***"
   fi

else

   echo "*** No compatible version of WRF physics tables found; attempting to download compatible tables ***"

fi

if [ ! -d physics_wrf/files ]; then
   mkdir -p physics_wrf/files
fi


#
# Try using 'git'
#
which git
if [ $? -eq 0 ]; then
   echo "*** Trying git to obtain WRF physics tables ***"
   git clone git://github.com/${github_org}/iAMAS-Data.git
   if [ $? -eq 0 ]; then
      cd iAMAS-Data
      git checkout v${iamas_vers}
      if [ $? -ne 0 ]; then
         echo "*** iAMAS version-specific tag not found; trying the master branch ***"
      else
         echo "*** Found v${iamas_vers} tag ***"
      fi
      cd ..
      mv iAMAS-Data/atmosphere/physics_wrf/files/* physics_wrf/files
      rm -rf iAMAS-Data

      check_compatibility
      if [ $? -eq 1 ]; then
         echo "*** Successfully obtained compatible versions of WRF physics tables ***"
         exit 0
      fi
   else
      echo "*** Failed to obtain WRF physics tables using git ***"
   fi
else
   echo "*** git not in path ***"
fi


#
# Try using 'svn'
#
which svn
if [ $? -eq 0 ]; then
   echo "*** Trying svn to obtain WRF physics tables ***"
   branch=v${iamas_vers}
   svn checkout --non-interactive --trust-server-cert https://github.com/${github_org}/iAMAS-Data.git/tags/${branch}
   if [ $? -ne 0 ]; then
      echo "*** iAMAS version-specific tag not found; trying the trunk ***"
      branch=trunk
      svn checkout --non-interactive --trust-server-cert https://github.com/${github_org}/iAMAS-Data.git/${branch}
   else
      echo "*** Found v${iamas_vers} tag ***"
   fi
   if [ $? -eq 0 ]; then
      mv ${branch}/atmosphere/physics_wrf/files/* physics_wrf/files
      rm -rf ${branch}
      check_compatibility
      if [ $? -eq 1 ]; then
         echo "*** Successfully obtained compatible versions of WRF physics tables ***"
         exit 0
      fi
   else
      echo "*** Failed to obtain WRF physics tables using svn ***"
   fi
else
   echo "*** svn not in path ***"
fi


#
# Try using 'curl'
#
which curl
if [ $? -eq 0 ]; then
   echo "*** Trying curl to obtain WRF physics tables ***"
   branch=${iamas_vers}
   curl -sf -o iAMAS-Data.tar.gz https://codeload.github.com/${github_org}/iAMAS-Data/tar.gz/v${branch}
   if [ $? -ne 0 ]; then
      echo "*** iAMAS version-specific tar file not found; trying the master tar file ***"
      branch=master
      curl -sf -o iAMAS-Data.tar.gz https://codeload.github.com/${github_org}/iAMAS-Data/tar.gz/${branch}
   else
      echo "*** Found v${iamas_vers} tar file ***"
   fi
   if [ $? -eq 0 ]; then
      which tar
      if [ $? -eq 0 ]; then
         tar -xzf iAMAS-Data.tar.gz
         mv iAMAS-Data-${branch}/atmosphere/physics_wrf/files/* physics_wrf/files
         rm -rf iAMAS-Data.tar.gz iAMAS-Data-${branch}

         check_compatibility
         if [ $? -eq 1 ]; then
            echo "*** Successfully obtained compatible versions of WRF physics tables ***"
            exit 0
         fi
      else
         echo "*** tar not in path -- unable to extract WRF physics tables ***"
         rm -rf iAMAS-Data.tar.gz
      fi
   else
      echo "*** Failed to obtain WRF physics tables using curl ***"
   fi
else
   echo "*** curl not in path ***"
fi


echo "***************************************************************"
echo "Unable to obtain WRF physics tables using git, svn, or curl."
echo "This may be because 'git', 'svn', and 'curl' are not installed,"
echo "or it could be due to network connectivity problems."
echo " "
echo "Please see src/core_atmosphere/physics/checkout_data_files.sh"
echo "for suggestions on how to remedy this issue."
echo "***************************************************************"

exit 1
