#!/bin/bash

function build_opm_grid {
  # Build ERT
  pushd .
  mkdir -p $WORKSPACE/deps/ert
  cd $WORKSPACE/deps/ert
  git init .
  git remote add origin https://github.com/Ensembles/ert
  git fetch --depth 1 origin $ERT_REVISION:branch_to_build
  test $? -eq 0 || exit 1
  git checkout branch_to_build
  popd

  pushd .
  mkdir -p serial/build-ert
  cd serial/build-ert
  cmake $WORKSPACE/deps/ert/devel -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$WORKSPACE/serial/install
  test $? -eq 0 || exit 1
  cmake --build . --target install
  test $? -eq 0 || exit 1
  popd

  # Build opm-common
  pushd .
  mkdir -p $WORKSPACE/deps/opm-common
  cd $WORKSPACE/deps/opm-common
  git init .
  git remote add origin https://github.com/OPM/opm-common
  git fetch --depth 1 origin $OPM_COMMON_REVISION:branch_to_build
  test $? -eq 0 || exit 1
  git checkout branch_to_build
  popd
  source $WORKSPACE/deps/opm-common/jenkins/build-opm-module.sh

  pushd .
  mkdir serial/build-opm-common
  cd serial/build-opm-common
  build_module "-DCMAKE_INSTALL_PREFIX=$WORKSPACE/serial/install" 0 $WORKSPACE/deps/opm-common
  test $? -eq 0 || exit 1
  popd

  # Build opm-parser
  clone_and_build_module opm-parser "-DCMAKE_PREFIX_PATH=$WORKSPACE/serial/install -DCMAKE_INSTALL_PREFIX=$WORKSPACE/serial/install" $OPM_PARSER_REVISION $WORKSPACE/serial
  test $? -eq 0 || exit 1

  # Build opm-material
  clone_and_build_module opm-material "-DCMAKE_PREFIX_PATH=$WORKSPACE/serial/install -DCMAKE_INSTALL_PREFIX=$WORKSPACE/serial/install" $OPM_MATERIAL_REVISION $WORKSPACE/serial
  test $? -eq 0 || exit 1

  # Build opm-core
  clone_and_build_module opm-core "-DCMAKE_PREFIX_PATH=$WORKSPACE/serial/install -DCMAKE_INSTALL_PREFIX=$WORKSPACE/serial/install" $OPM_CORE_REVISION $WORKSPACE/serial
  test $? -eq 0 || exit 1

  # Build opm-grid
  pushd .
  mkdir serial/build-opm-grid
  cd serial/build-opm-grid
  build_module "-DCMAKE_PREFIX_PATH=$WORKSPACE/serial/install" 1 $WORKSPACE
  test $? -eq 0 || exit 1
  popd
}
