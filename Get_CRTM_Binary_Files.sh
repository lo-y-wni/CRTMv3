#https://bin.ssec.wisc.edu/pub/s4/CRTM/fix_REL-3.1.1.2.tgz  (use this for jedi and stand-alone, some files have changed).

# This script is used to manually download the tarball of binary and netcdf coefficient files.
# The same files also download automatically during the cmake step, so you don't have to actually run this manually. 

foldername="fix_REL-3.1.1.2"
checksum=58e0a5c698a438a31dc4914fcda39846
filename="${foldername}.tgz"
echo "$filename"

download_url=https://bin.ssec.wisc.edu/pub/s4/CRTM/$filename
md5sum_url="${download_url}.md5sum"


# Check if the fix directory already exists, indicating that no download is needed.
if [ -d "fix/" ]; then #fix directory exists
    echo "fix/ already exists, doing nothing."
    exit 0
fi

# If var $CRTM_BINARY_FILES_TARBALL is set, confirm the checksum then
# update the "filename" var.
if [ -f "$CRTM_BINARY_FILES_TARBALL" ]; then
    local_checksum=$(md5sum $CRTM_BINARY_FILES_TARBALL | cut -d ' ' -f 1)
    if [ "${local_checksum}" = "${checksum}" ]; then
        filename=$CRTM_BINARY_FILES_TARBALL
    else
        # Exit with failure.
        echo "Var CRTM_BINARY_FILES_TARBALL is set but does not match expected"
        echo "checksum of $checksum. Confirm you have the correct tarball and"
        echo "try again, or unset this variable"
        exit 1
    fi
fi

# If the file is not present in the pwd (or at the location set by
# $CRTM_BINARY_FILES_TARBALL) download the coefficients file.
if ! test -f "$filename"; then
    echo "Downloading $filename, please wait about 7 minutes (7 GB tar file: sorry!)"
    wget $download_url # CRTM binary files, add "-q" to suppress output.
fi

# Extract the file to pwd and write the checksum.
untar the file and move directory to fix
tar -zxvf $filename -C "${PWD}"
mkdir -p fix
mv $foldername/fix/* fix/.
rm -rf $foldername

echo "Completed."
