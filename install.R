local_lib_path = "~/R_libs/"

install.packages('minpack.lm', lib=local_lib_path)
install.packages('argparse', lib=local_lib_path)
install.packages('.', repos=NULL, type="source", lib=local_lib_path)
