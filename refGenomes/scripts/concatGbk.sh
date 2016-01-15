for dir in */
do
    cd $dir
    echo $dir
    cat *.gbk > ${dir%/}.gbk
    mv ${dir%/}.gbk ../${dir%/}.gbk
    cd ../
done
