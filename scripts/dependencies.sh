#!/bin/bash
set -e

export LD_LIBRARY_PATH=.:`cat /etc/ld.so.conf.d/* | grep -v -E "#" | tr "\\n" ":" | sed -e "s/:$//g"`
if [[ $(grep -cF trusty /etc/lsb-release) > 0 ]]; then
  sudo add-apt-repository -y ppa:kevinkreiser/libsodium
  sudo add-apt-repository -y ppa:kevinkreiser/libpgm
  sudo add-apt-repository -y ppa:kevinkreiser/zeromq3
  sudo add-apt-repository -y ppa:kevinkreiser/czmq
fi
sudo add-apt-repository -y ppa:kevinkreiser/prime-server
sudo apt-get update
sudo apt-get install -y autoconf automake pkg-config libtool make pkg-config gcc g++ lcov libboost1.54-all-dev libprime-server0.6.0-dev

#clone async
mkdir -p deps
for dep in midgard baldr sif; do
	git clone --depth=1 --recurse --single-branch https://github.com/valhalla/$dep.git deps/$dep &
done
wait

#build sync
for dep in midgard baldr sif; do
	pushd deps/$dep
	./autogen.sh
	./configure CPPFLAGS="-DBOOST_SPIRIT_THREADSAFE -DBOOST_NO_CXX11_SCOPED_ENUMS"
	make -j$(nproc)
	sudo make install
	popd
done
