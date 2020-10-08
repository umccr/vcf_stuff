# `--single-version-externally-managed --root=/` is added to make setuptools avoid downloading dependencies. This resolves the conda-build
#       error "Setuptools downloading is disabled in conda build. Be sure to add all dependencies in the meta.yaml"
$PYTHON setup.py install --single-version-externally-managed --root=/
chmod -R o+r $PREFIX/lib/python*/site-packages/*
