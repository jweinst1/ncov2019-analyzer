echo "Running the covid-19 shell script"
mkdir cov19_build
cd cov19_build
cmake ..
make && make test
cd ..
rm -rf cov19-build
echo "Shell script complete"