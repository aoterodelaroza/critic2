set(TESTS
  001_crystalenv-simple
  002_crystalenv-arithm
  003_crystalenv-opts
  004_crystalenv-tricky1
  005_crystalenv-tricky2
  006_crystalenv-tricky3
  007_crystal_abinit
  008_crystal_crystal
  009_crystal_elk
  010_crystal_cifallspg
  011_crystal_cifone
  012_crystal_cifallblock
  013_crystal_cifoneblock
  )

runtests(${TESTS})
