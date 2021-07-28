rm -rf cressp/data/pdb/makeblastdb_out
rm -f /home/merit_an/git/cressp/cressp/data/pdb/rcsb_pdb.fa
rm -f /home/merit_an/git/cressp/cressp/data/pdb/swiss_model.fa
git add .
git add -A
git commit -am "update"
git push
