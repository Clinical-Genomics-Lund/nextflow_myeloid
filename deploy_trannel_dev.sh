DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

PIPELINE_DEST="/trannel/proj/varlociraptor/"
#DEST_HOST="rs-fs1.lunarc.lu.se"


# Copy pipeline script
cp $DIR/main.nf $PIPELINE_DEST
cp $DIR/varli_viktor.nf $PIPELINE_DEST


# Copy configuration file
cp $DIR/configs/nextflow.trannel.config $PIPELINE_DEST/nextflow.config

# Copy other files
cp -r $DIR/bin $PIPELINE_DEST

git rev-parse HEAD > git.hash
cp $DIR/git.hash $PIPELINE_DEST