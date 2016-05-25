for folder in *;
do mv "$folder" "${folder/KBase_joshamilton:1464097239222_/}";
done

for folder in *;
do mv "$folder" "${folder/\.model_KBaseFBA.FBAModel-8.0/}";
done

for file in */*;
do mv "$file" "${file/joshamilton:1464097239222-/}";
done

for file in */*;
do mv "$file" "${file/joshamilton:1464097239222_/}";
done

for file in */*;
do mv "$file" "${file/.model-SBML/}";
done

for file in */*;
do mv "$file" "${file/.model_FBAModel/}";
done
